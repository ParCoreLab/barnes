/*************************************************************************/
/*                                                                       */
/*  Copyright (c) 1994 Stanford University                               */
/*                                                                       */
/*  All rights reserved.                                                 */
/*                                                                       */
/*  Permission is given to use, copy, and modify this software for any   */
/*  non-commercial purpose as long as this copyright notice is not       */
/*  removed.  All other uses, including redistribution in whole or in    */
/*  part, are forbidden without prior written permission.                */
/*                                                                       */
/*  This software is provided with absolutely no warranty and no         */
/*  support.                                                             */
/*                                                                       */
/*************************************************************************/
/*
Usage: BARNES <options> < inputfile

Command line options:

    -h : Print out input file description

    Input parameters should be placed in a file and redirected through
    standard input.  There are a total of twelve parameters, and all of
    them have default values.

    1) infile (char*) : The name of an input file that contains particle
       data.

       The format of the file is:
         a) An int representing the number of particles in the distribution
         b) An int representing the dimensionality of the problem (3-D)
         c) A double representing the current time of the simulation
         d) Doubles representing the masses of all the particles
         e) A vector (length equal to the dimensionality) of doubles
            representing the positions of all the particles
         f) A vector (length equal to the dimensionality) of doubles
            representing the velocities of all the particles

       Each of these numbers can be separated by any amount of whitespace.
    2) nbody (int) : If no input file is specified (the first line is
       blank), this number specifies the number of particles to generate
       under a plummer model.  Default is 16384.
    3) seed (int) : The seed used by the random number generator.
       Default is 123.
    4) outfile (char*) : The name of the file that snapshots will be
       printed to. This feature has been disabled in the SPLASH release.
       Default is NULL.
    5) dtime (double) : The integration time-step.
       Default is 0.025.
    6) eps (double) : The usual potential softening
       Default is 0.05.
    7) tol (double) : The cell subdivision tolerance.
       Default is 1.0.
    8) fcells (double) : Number of cells created = fcells * number of
       leaves.
       Default is 2.0.
    9) fleaves (double) : Number of leaves created = fleaves * nbody.
       Default is 0.5.
    10) tstop (double) : The time to stop integration.
       Default is 0.075.
    11) dtout (double) : The data-output interval.
       Default is 0.25.
    12) NPROC (int) : The number of processors.
       Default is 1.
*/

MAIN_ENV

#include <chrono>
#include <cassert>
#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <bits/stdc++.h>
#include <unordered_map>
#include <fcntl.h>

#define global  /* nada */

#include "stdinc.h"
#include "cha.h"
#include "topology.h"
#include "load.h"
#include "grav.h"
#include "constants.h"
#define NTHREADS 28
//#include "hash.h"
#define CACHE_LINE_SZ (64)
#define ALIGN_TO_CACHE_LINE(addr) ((uint64_t)(addr) & (~(CACHE_LINE_SZ-1)))
#define OFFSET_OF_CACHE_LINE(addr) ((uint64_t)(addr) >> 6)
#define HASHTABLE_SIZE (100000)
#define SAMPLING_PERIOD (100000)
// begin
//#if 0
//using namespace std;

class HashTable
{
        int INDEX_NUM;    // No. of buckets

        // Pointer to an array containing buckets
        std::pair<std::atomic<uint32_t>,std::list<std::tuple<long unsigned int, int, std::pair<int, int>>>> *table;
        public:
        HashTable(int index_num);  // Constructor
        ~HashTable();  // Destructor

        // inserts a key into hash table
        void insertItem(long unsigned int index, long unsigned int cacheLine, int cha);

        bool findItem(long unsigned int index, long unsigned int cacheLine);

        void wait(long unsigned int index);

        void signal(long unsigned int index);

        std::tuple<long unsigned int, int, std::pair<int, int>>& getItem(long unsigned int index, long unsigned int cacheLine);

        // hash function to map values to key
        long unsigned int hashFunction(long unsigned int cacheLine) {
                return (cacheLine % INDEX_NUM);
        }

        void displayHash();
};

HashTable::HashTable(int index_num)
{
        this->INDEX_NUM = index_num;
        table = new std::pair<std::atomic<uint32_t>,std::list<std::tuple<long unsigned int, int, std::pair<int, int>>>>[INDEX_NUM];
        for(int i = 0; i < INDEX_NUM; i++)
                table[i].first = 1;
}

HashTable::~HashTable()
{
        delete[] table;
}

void HashTable::insertItem(long unsigned int index, long unsigned int cacheLine, int cha)
{
        //long unsigned int index = hashFunction(cacheLine);
        table[index].second.push_back({cacheLine, cha, {0, 0}});
}

bool HashTable::findItem(long unsigned int index, long unsigned int cacheLine)
{
        // get the hash index of key
        //long unsigned int index = hashFunction(cacheLine);

        auto it = table[index].second.begin();
        while (it != table[index].second.end())
        {
                if (std::get<0>(*it) == cacheLine)
                        break;
        }

        // if key is found in hash table, remove it
        if (it != table[index].second.end())
                return true;
        return false;
}
//#endif

void HashTable::wait(long unsigned int index)
{
        //while(table[index].first <= 0);
        //table[index].first--;
        auto oldval = table[index].first.load();
        while (oldval == 0 || !table[index].first.compare_exchange_strong(oldval, oldval - 1)) {
           oldval = table[index].first.load();
           //std::cerr << "in wait, val: " << oldval << "\n";
        }
}

void HashTable::signal(long unsigned int index)
{
        auto val = table[index].first.fetch_add(1, std::memory_order_seq_cst);
        //std::cerr << "in signal, val: " << val << "\n";
}

std::tuple<long unsigned int, int, std::pair<int, int>>& HashTable::getItem(long unsigned int index, long unsigned int cacheLine)
{
        // get the hash index of key
        //long unsigned int index = hashFunction(cacheLine);

        // find the key in (index)th list
        auto it = table[index].second.begin();
        while (it != table[index].second.end())
        {
                if (std::get<0>(*it) == cacheLine)
                        break;
        }

        // if key is found in hash table, remove it
        if (it != table[index].second.end()) {
                return *it;//table[index];
        }
        std::tuple<long unsigned int, int, std::pair<int, int>> empty_tuple({0, 0, {0, 0}});
        return empty_tuple;
}

// function to display hash table
void HashTable::displayHash() {
        for (int i = 0; i < INDEX_NUM; i++) {
                std::cout << i;
                for (auto x : table[i].second)
                        //cout << " --> cache line: " << x.first << ", a: " << x.second.first << ", b: " << x.second.second;
                        std::cout << " --> cache line: " << std::get<0>(x) << ", cha: " << std::get<1>(x) << ", a: " << std::get<2>(x).first << ", b: " << std::get<2>(x).second;
                std::cout << std::endl;
        }
}
//#endif
// end

std::map<long, std::multiset<double *>> threadid_addresses_map;
HashTable comm_map(HASHTABLE_SIZE);
std::vector<std::vector<std::pair<int, std::unordered_map<int, int>>>> comm_matrix(NTHREADS, std::vector<std::pair<int, std::unordered_map<int, int>>> (NTHREADS, std::pair<int, std::unordered_map<int, int>> {}));

pthread_spinlock_t map_spinlock;

string defv[] = {                 /* DEFAULT PARAMETER VALUES              */
    /* file names for input/output                                         */
    "in=",                        /* snapshot of initial conditions        */
    "out=",                       /* stream of output snapshots            */

    /* params, used if no input specified, to make a Plummer Model         */
    "nbody=16384",                /* number of particles to generate       */
    "seed=123",                   /* random number generator seed          */

    /* params to control N-body integration                                */
    "dtime=0.025",                /* integration time-step                 */
    "eps=0.05",                   /* usual potential softening             */
    "tol=1.0",                    /* cell subdivision tolerence            */
    "fcells=2.0",                 /* cell allocation parameter             */
    "fleaves=0.5",                 /* leaf allocation parameter             */

    "tstop=0.075",                 /* time to stop integration              */
    "dtout=0.25",                 /* data-output interval                  */

    "NPROC=1",                    /* number of processors                  */
};

/* The more complicated 3D case */
#define NUM_SOCKETS 2
#define NUM_CHA_BOXES 28
#define NUM_DIRECTIONS 32
#define BRC_FUC 0
#define BRC_FRA 1
#define BRA_FDA 2
#define BRA_FRC 3
#define BLC_FDC 4
#define BLC_FLA 5
#define BLA_FUA 6
#define BLA_FLC 7
#define BUC_FUA 8
#define BUC_FLC 9
#define BUA_FUC 10
#define BUA_FRA 11
#define BDC_FDA 12
#define BDC_FRC 13
#define BDA_FDC 14
#define BDA_FLA 15

#define FRC_BUC 16
#define FRC_BRA 17
#define FRA_BDA 18
#define FRA_BRC 19
#define FLC_BDC 20
#define FLC_BLA 21
#define FLA_BUA 22
#define FLA_BLC 23
#define FUC_BUA 24
#define FUC_BLC 25
#define FUA_BUC 26
#define FUA_BRA 27
#define FDC_BDA 28
#define FDC_BRC 29
#define FDA_BDC 30
#define FDA_BLA 31

#define CACHELINE_SIZE 64

unsigned short rdtsc() {
        unsigned short c;
        __asm__("rdtsc\n" : "=a" (c));
        return c;
}

thread_local uint64_t sampling_counter = 0;
thread_local uint64_t next_sampling_iteration = 0;

//int findCha(const double* val);
int findCha(uint64_t val);

void inc_comm(int tid,
                //HashTable& comm_map,
                //std::vector<std::vector<std::pair<int, std::unordered_map<int, int>>>>& comm_matrix,
                //const int& cha,
                const long unsigned int& addr)
{
//#if 0 
        if(sampling_counter++ < next_sampling_iteration) {
                //cerr << "inc_comm discarded\n";
                return;
        }
//#endif
        next_sampling_iteration = SAMPLING_PERIOD - 256 + rdtsc() % 256 + next_sampling_iteration;
//#endif
#if 0
        if(tid == 0)
                cerr << "next_sampling_iteration: " << next_sampling_iteration << endl;
#endif  
        uint64_t line = OFFSET_OF_CACHE_LINE(addr);
        uint64_t hash_idx = comm_map.hashFunction(line);
        int sh = 1;
        int a;
        int b;
	int cha;
        
        comm_map.wait(hash_idx);
        if(!comm_map.findItem(hash_idx, line))  {
		cha = findCha(line);
                comm_map.insertItem(hash_idx, line, cha);
        }
        std::tuple<long unsigned int, int, std::pair<int, int>>& comm_map_element = comm_map.getItem(hash_idx, line);
	cha = a = std::get<1>(comm_map_element);
        //a = comm_map_element.second.first;
	a = std::get<2>(comm_map_element).first;   
        //b = comm_map_element.second.second;
        b = std::get<2>(comm_map_element).second;
        
        if (a == 0 && b == 0)
                sh = 0;
        if (a != 0 && b != 0)
                sh = 2;
        switch (sh) {
                case 0: // no one accessed line before, store accessing thread in pos 0
                        //comm_map[line].first = tid+1;
                        //comm_map[line].first = tid+1;
                        //comm_map_element.second.first = tid+1;
                        std::get<2>(comm_map_element).first = tid+1;
                        //mutex_map[hash_idx].unlock();
                        comm_map.signal(hash_idx);
                        //mtx->unlock();
                        //global_mtx.unlock();
                        //comm_map[line].first = tid+1;
                        break;

                case 1: // one previous access => needs to be in pos 0
                        // if (a != tid+1) {
                        //inc_comm(tid, a);
                        //if (tid!=a-1) {
                        //      if(tid < a-1) {
                        //comm_map[line].first = tid+1;
                        //comm_map_element.second.first = tid+1;
                        std::get<2>(comm_map_element).first = tid+1;
                        //comm_map[line].second = a;
                        //comm_map_element.second.second = a;
                        std::get<2>(comm_map_element).second = a;
                        //mutex_map[hash_idx].unlock();
                        //comm_map.signal(hash_idx);
                        //mtx->unlock();
                        //global_mtx.unlock();
                        comm_matrix[tid][a-1].first++;
                        comm_matrix[tid][a-1].second[cha]++;
			comm_map.signal(hash_idx);
                        //      } else {
                        //              comm_matrix[a-1][tid].first++;
                        //              comm_matrix[a-1][tid].second[cha]++;
                        //      }
                        //}
                        //#pragma omp critical
                        //{
                        //comm_map[line].first = tid+1;
                        //comm_map[line].second = a;
                        //}
                        // }
                        break;
                        
                case 2: // two previous accesses
                        // if (a != tid+1 && b != tid+1) {
                        //comm_map[line].first = tid+1;
                        //comm_map_element.second.first = tid+1;
                        std::get<2>(comm_map_element).first = tid+1;
                        //comm_map[line].second = a;
                        //comm_map_element.second.second = a;
                        std::get<2>(comm_map_element).second = a;
                        //mutex_map[hash_idx].unlock();
                        //comm_map.signal(hash_idx);
                        //mtx->unlock();
                        //global_mtx.unlock();
                        if (tid!=a-1) {
                                //if(tid < a-1) {
                                comm_matrix[tid][a-1].first++;
                                comm_matrix[tid][a-1].second[cha]++;
                                //} else {
                                //      comm_matrix[a-1][tid].first++;
                                //      comm_matrix[a-1][tid].second[cha]++; 
                                //}
                        }
                        if (tid!=b-1) {
                                //if(tid < b-1) {
                                comm_matrix[tid][b-1].first++;
                                comm_matrix[tid][b-1].second[cha]++;
                                //} else {
                                //      comm_matrix[b-1][tid].first++;
                                //      comm_matrix[b-1][tid].second[cha]++;
                                //}
                        }
			comm_map.signal(hash_idx);
                        //#pragma omp critical
                        //{
                        //comm_map[line].first = tid+1;
                        //comm_map[line].second = a;
                        //}
                        // } else if (a == tid+1) {
                        //  inc_comm(tid, b);
                        // } else if (b == tid+1) {
                        //  inc_comm(tid, a);
                        //  commmap[line].first = tid+1;
                        //  commmap[line].second = a;
                        // }

                        break;
        }
        //mtx->unlock();                        
        //#endif
        //__sync_synchronize();
        //#endif
        //std::cerr << "This part is executed 1\n"; 
}

//static const long CHA_MSR_PMON_CTRL_BASE = 0x0E01L;

enum class TrafficType
{
	Data,
	Address,
	Acknowledge,
	Invalidate
};

std::string enumToStr(TrafficType traffic_type)
{
	std::string res;
	switch (traffic_type)
	{
		case TrafficType::Data:
			res = "data";
			break;
		case TrafficType::Address:
			res = "address";
			break;
		case TrafficType::Acknowledge:
			res = "acknowledge";
			break;
		case TrafficType::Invalidate:
			res = "invalidate";
			break;
		default:
			res = "unknown";
			break;
	}

	return res;
}

std::map<int, int> getMsrFds()
{
        // SPDLOG_TRACE(__PRETTY_FUNCTION__);

        std::map<int, int> fd_map;

        char filename[100];

        auto logical_core_count = getCoreCount();
        // SPDLOG_TRACE("logical core count: {}", logical_core_count);

        for(auto i = 0; i < logical_core_count; ++i) {
            sprintf(filename, "/dev/cpu/%d/msr",i);

            int fd = open(filename, O_RDWR);

            if(fd >= 0) {
                // SPDLOG_TRACE("Opened fd: {}", fd);
                fd_map.insert({i, fd});
            } else if(fd == -1) {
                // SPDLOG_ERROR("error on open(): {}", strerror(errno));
            }
        }

        for(const auto& p : fd_map) {
            // SPDLOG_TRACE("MSR fd of core {}: {}", p.first, p.second);
        }

        return fd_map;
}

void setAllUncoreRegisters(const std::vector<unsigned int>& vals)
{
        // SPDLOG_TRACE(__PRETTY_FUNCTION__);

        // SPDLOG_TRACE("Setting all uncore registers with values: ");
        for(auto val : vals) {
            // SPDLOG_TRACE("{:x}", val);
        }

        int processor_in_socket[NUM_SOCKETS];
        int logical_core_count = getCoreCount();
        processor_in_socket[0] = 0;
        processor_in_socket[1] = logical_core_count - 1;

        auto msr_fds = getMsrFds();

        for(auto socket = 0; socket < NUM_SOCKETS; ++socket) {
            for(auto cha = 0; cha < NUM_CHA_BOXES; ++cha) {
                int core = processor_in_socket[socket];

                for(auto i = 0u; i < vals.size(); ++i) {
                    uint64_t val = vals[i];
                    uint64_t offset = CHA_MSR_PMON_CTRL_BASE + (0x10 * cha) + i;

                    ssize_t rc64 = pwrite(msr_fds[core], &val, sizeof(val), offset);
                    if(rc64 == 8) {
                        // SPDLOG_TRACE("Configuring socket {}, CHA {}, by writing 0x{:x} to core {} (fd: {}), offset 0x{:x}.",
//                                      socket, cha, val, core, msr_fds[core], offset);
                    } else {
                        // SPDLOG_ERROR("Error writing all data to MSR device on core {}, written {} bytes.", core, rc64);
                    }
                }
            }
        }

        /// it is important to close this as well, otherwise we will have fd leak.
        // SPDLOG_TRACE("closing file descriptors of MSRs.");
        for(const auto& p : msr_fds) {
            int cpu = p.first;
            int to_be_closed = p.second;
            // SPDLOG_TRACE("closing fd {} of cpu {}.", to_be_closed, cpu);
            ::close(to_be_closed);
            // SPDLOG_TRACE("closed fd: {}", to_be_closed);
        }
}

std::vector<std::vector<uint64_t>> storeTraffic(TrafficType traffic_type)
{
	std::vector<std::vector<uint64_t>> res;

	auto msr_fds = getMsrFds();

	unsigned int left_val = 0;
	unsigned int right_val = 0;
	unsigned int up_val = 0;
	unsigned int down_val = 0;

	if (traffic_type == TrafficType::Data)
	{
		left_val = LEFT_BL_READ;
		right_val = RIGHT_BL_READ;
		up_val = UP_BL_READ;
		down_val = DOWN_BL_READ;
	}
	else if (traffic_type == TrafficType::Address)
	{
		left_val = LEFT_AD_READ;
		right_val = RIGHT_AD_READ;
		up_val = UP_AD_READ;
		down_val = DOWN_AD_READ;
	}
	else if (traffic_type == TrafficType::Acknowledge)
	{
		left_val = LEFT_AK_READ;
		right_val = RIGHT_AK_READ;
		up_val = UP_AK_READ;
		down_val = DOWN_AK_READ;
	}
	else if (traffic_type == TrafficType::Invalidate)
	{
		left_val = LEFT_IV_READ;
		right_val = RIGHT_IV_READ;
		up_val = UP_IV_READ;
		down_val = DOWN_IV_READ;
	}
	else
	{
		std::cerr << "unexpected traffic type!";
		std::exit(1);
	}

	std::vector<unsigned int> vals = {left_val, right_val, up_val, down_val, FILTER0_OFF, FILTER1_OFF};
	setAllUncoreRegisters(vals);

	/// popping up filter values since we will not read from them.
	vals.pop_back();
	vals.pop_back();

	//SPDLOG_DEBUG("---------------- FIRST READINGS ----------------");
	for (int socket = 0; socket < /*NUM_SOCKETS*/ 1; ++socket)
	{                  // we will work on socket-0
		long core = 0; // this line will only work for socket-0!

		for (int cha = 0; cha < NUM_CHA_BOXES; ++cha)
		{
			std::vector<uint64_t> direction_vals;

			for (int i = 0; i < vals.size(); ++i)
			{
				uint64_t msr_val;
				uint64_t msr_num = CHA_MSR_PMON_CTR_BASE + (CHA_BASE * cha) + i;
				//SPDLOG_DEBUG("Executing pread() --> fd: {}, offset: {:x}", msr_fds[core], msr_num);
				ssize_t rc64 = pread(msr_fds[core], &msr_val, sizeof(msr_val), msr_num);
				if (rc64 != sizeof(msr_val))
				{
					//SPDLOG_ERROR("EXIT FAILURE. rc64: {}", rc64);
					//SPDLOG_ERROR("error: {}", strerror(errno));
					exit(EXIT_FAILURE);
				}
				else
				{
					direction_vals.push_back(msr_val);
					//SPDLOG_DEBUG("Read {} from socket {}, CHA {} on core {}, offset 0x{:x}.", msr_val, socket, cha, core, msr_num);
				}
			}

			res.push_back(direction_vals);
		}
	}

	//SPDLOG_DEBUG("closing file descriptors of MSRs.");
	for (const auto &p : msr_fds)
	{
		int cpu = p.first;
		int to_be_closed = p.second;
		//SPDLOG_DEBUG("closing fd {} of cpu {}.", to_be_closed, cpu);
		::close(to_be_closed);
	}

	return res;
}

enum class BenchmarkOption
{
	None,
	OnlyDefault,
	OnlyChaAware,
	Both
};


static long Child_Sequence[NUM_DIRECTIONS][NSUB] =
{
	{ 2, 5, 6, 1, 0, 3, 4, 7},  /* BRC_FUC */
 { 2, 5, 6, 1, 0, 7, 4, 3},  /* BRC_FRA */
 { 1, 6, 5, 2, 3, 0, 7, 4},  /* BRA_FDA */
 { 1, 6, 5, 2, 3, 4, 7, 0},  /* BRA_FRC */
 { 6, 1, 2, 5, 4, 7, 0, 3},  /* BLC_FDC */
 { 6, 1, 2, 5, 4, 3, 0, 7},  /* BLC_FLA */
 { 5, 2, 1, 6, 7, 4, 3, 0},  /* BLA_FUA */
 { 5, 2, 1, 6, 7, 0, 3, 4},  /* BLA_FLC */
 { 1, 2, 5, 6, 7, 4, 3, 0},  /* BUC_FUA */
 { 1, 2, 5, 6, 7, 0, 3, 4},  /* BUC_FLC */
 { 6, 5, 2, 1, 0, 3, 4, 7},  /* BUA_FUC */
 { 6, 5, 2, 1, 0, 7, 4, 3},  /* BUA_FRA */
 { 5, 6, 1, 2, 3, 0, 7, 4},  /* BDC_FDA */
 { 5, 6, 1, 2, 3, 4, 7, 0},  /* BDC_FRC */
 { 2, 1, 6, 5, 4, 7, 0, 3},  /* BDA_FDC */
 { 2, 1, 6, 5, 4, 3, 0, 7},  /* BDA_FLA */

 { 3, 4, 7, 0, 1, 2, 5, 6},  /* FRC_BUC */
 { 3, 4, 7, 0, 1, 6, 5, 2},  /* FRC_BRA */
 { 0, 7, 4, 3, 2, 1, 6, 5},  /* FRA_BDA */
 { 0, 7, 4, 3, 2, 5, 6, 1},  /* FRA_BRC */
 { 7, 0, 3, 4, 5, 6, 1, 2},  /* FLC_BDC */
 { 7, 0, 3, 4, 5, 2, 1, 6},  /* FLC_BLA */
 { 4, 3, 0, 7, 6, 5, 2, 1},  /* FLA_BUA */
 { 4, 3, 0, 7, 6, 1, 2, 5},  /* FLA_BLC */
 { 0, 3, 4, 7, 6, 5, 2, 1},  /* FUC_BUA */
 { 0, 3, 4, 7, 6, 1, 2, 5},  /* FUC_BLC */
 { 7, 4, 3, 0, 1, 2, 5, 6},  /* FUA_BUC */
 { 7, 4, 3, 0, 1, 6, 5, 2},  /* FUA_BRA */
 { 4, 7, 0, 3, 2, 1, 6, 5},  /* FDC_BDA */
 { 4, 7, 0, 3, 2, 5, 6, 1},  /* FDC_BRC */
 { 3, 0, 7, 4, 5, 6, 1, 2},  /* FDA_BDC */
 { 3, 0, 7, 4, 5, 2, 1, 6},  /* FDA_BLA */
};

static long Direction_Sequence[NUM_DIRECTIONS][NSUB] =
{
	{ FRC_BUC, BRA_FRC, FDA_BDC, BLA_FUA, BUC_FLC, FUA_BUC, BRA_FRC, FDA_BLA },
 /* BRC_FUC */
 { FRC_BUC, BRA_FRC, FDA_BDC, BLA_FUA, BRA_FDA, FRC_BRA, BUC_FUA, FLC_BDC },
 /* BRC_FRA */
 { FRA_BDA, BRC_FRA, FUC_BUA, BLC_FDC, BDA_FLA, FDC_BDA, BRC_FRA, FUC_BLC },
 /* BRA_FDA */
 { FRA_BDA, BRC_FRA, FUC_BUA, BLC_FDC, BUC_FLC, FUA_BUC, BRA_FRC, FDA_BLA },
 /* BRA_FRC */
 { FLC_BDC, BLA_FLC, FUA_BUC, BRA_FDA, BDC_FRC, FDA_BDC, BLA_FLC, FUA_BRA },
 /* BLC_FDC */
 { FLC_BDC, BLA_FLC, FUA_BUC, BRA_FDA, BLA_FUA, FLC_BLA, BDC_FDA, FRC_BUC },
 /* BLC_FLA */
 { FLA_BUA, BLC_FLA, FDC_BDA, BRC_FUC, BUA_FRA, FUC_BUA, BLC_FLA, FDC_BRC },
 /* BLA_FUA */
 { FLA_BUA, BLC_FLA, FDC_BDA, BRC_FUC, BLC_FDC, FLA_BLC, BUA_FUC, FRA_BDA },
 /* BLA_FLC */
 { FUC_BLC, BUA_FUC, FRA_BRC, BDA_FLA, BUA_FRA, FUC_BUA, BLC_FLA, FDC_BRC },
 /* BUC_FUA */
 { FUC_BLC, BUA_FUC, FRA_BRC, BDA_FLA, BLC_FDC, FLA_BLC, BUA_FUC, FRA_BDA },
 /* BUC_FLC */
 { FUA_BRA, BUC_FUA, FLC_BLA, BDC_FRC, BUC_FLC, FUA_BUC, BRA_FRC, FDA_BLA },
 /* BUA_FUC */
 { FUA_BRA, BUC_FUA, FLC_BLA, BDC_FRC, BRA_FDA, FRC_BRA, BUC_FUA, FLC_BDC },
 /* BUA_FRA */
 { FDC_BRC, BDA_FDC, FLA_BLC, BUA_FRA, BDA_FLA, FDC_BDA, BRC_FRA, FUC_BLC },
 /* BDC_FDA */
 { FDC_BRC, BDA_FDC, FLA_BLC, BUA_FRA, BUC_FLC, FUA_BUC, BRA_FRC, FDA_BLA },
 /* BDC_FRC */
 { FDA_BLA, BDC_FDA, FRC_BRA, BUC_FLC, BDC_FRC, FDA_BDC, BLA_FLC, FUA_BRA },
 /* BDA_FDC */
 { FDA_BLA, BDC_FDA, FRC_BRA, BUC_FLC, BLA_FUA, FLC_BLA, BDC_FDA, FRC_BUC },
 /* BDA_FLA */

 { BUC_FLC, FUA_BUC, BRA_FRC, FDA_BLA, FUC_BLC, BUA_FUC, FRA_BRC, BDA_FLA },
 /* FRC_BUC */
 { BUC_FLC, FUA_BUC, BRA_FRC, FDA_BLA, FRA_BDA, BRC_FRA, FUC_BUA, BLC_FDC },
 /* FRC_BRA */
 { BRA_FDA, FRC_BRA, BUC_FUA, FLC_BDC, FDA_BLA, BDC_FDA, FRC_BRA, BUC_FLC },
 /* FRA_BDA */
 { BRA_FDA, FRC_BRA, BUC_FUA, FLC_BDC, FRC_BUC, BRA_FRC, FDA_BDC, BLA_FUA },
 /* FRA_BRC */
 { BLC_FDC, FLA_BLC, BUA_FUC, FRA_BDA, FDC_BRC, BDA_FDC, FLA_BLC, BUA_FRA },
 /* FLC_BDC */
 { BLC_FDC, FLA_BLC, BUA_FUC, FRA_BDA, FLA_BUA, BLC_FLA, FDC_BDA, BRC_FUC },
 /* FLC_BLA */
 { BLA_FUA, FLC_BLA, BDC_FDA, FRC_BUC, FUA_BRA, BUC_FUA, FLC_BLA, BDC_FRC },
 /* FLA_BUA */
 { BLA_FUA, FLC_BLA, BDC_FDA, FRC_BUC, FLC_BDC, BLA_FLC, FUA_BUC, BRA_FDA },
 /* FLA_BLC */
 { BUC_FLC, FUA_BUC, BRA_FRC, FDA_BLA, FUA_BRA, BUC_FUA, FLC_BLA, BDC_FRC },
 /* FUC_BUA */
 { BUC_FLC, FUA_BUC, BRA_FRC, FDA_BLA, FLC_BDC, BLA_FLC, FUA_BUC, BRA_FDA },
 /* FUC_BLC */
 { BUA_FRA, FUC_BUA, BLC_FLA, FDC_BRC, FUC_BLC, BUA_FUC, FRA_BRC, BDA_FLA },
 /* FUA_BUC */
 { BUA_FRA, FUC_BUA, BLC_FLA, FDC_BRC, FRA_BDA, BRC_FRA, FUC_BUA, BLC_FDC },
 /* FUA_BRA */
 { BDC_FRC, FDA_BDC, BLA_FLC, FUA_BRA, FDA_BLA, BDC_FDA, FRC_BRA, BUC_FLC },
 /* FDC_BDA */
 { BDC_FRC, FDA_BDC, BLA_FLC, FUA_BRA, FRC_BUC, BRA_FRC, FDA_BDC, BLA_FUA },
 /* FDC_BRC */
 { BDA_FLA, FDC_BDA, BRC_FRA, FUC_BLC, FDC_BRC, BDA_FDC, FLA_BLC, BUA_FRA },
 /* FDA_BDC */
 { BDA_FLA, FDC_BDA, BRC_FRA, FUC_BLC, FLA_BUA, BLC_FLA, FDC_BDA, BRC_FUC },
 /* FDA_BLA */
};

int findCha(/*const double**/uint64_t val)
{
  // this part is changed wrt fluidanimate.
    return findCHAByHashing(reinterpret_cast<uintptr_t>(val)); // AYDIN: this is not &val, right?
}

int getMostAccessedCHA(int tid1,
                       int tid2,
                       std::multiset<std::tuple<int, int, int, int>, std::greater<std::tuple<int, int, int, int>>> ranked_cha_access_count_per_pair,
                       Topology topo)
{
    int max = 0;
    std::vector<int> considered_chas;
    std::map<int, bool> considered_chas_flag;

    auto it = ranked_cha_access_count_per_pair.begin();
    while (it != ranked_cha_access_count_per_pair.end())
    {
        // std::pair<int, int> tid_pair(std::get<2>(*it), std::get<3>(*it));
        if ((std::get<2>(*it) == tid1 && std::get<3>(*it) == tid2) || (std::get<3>(*it) == tid1 && std::get<2>(*it) == tid2))
        {
            // SPDLOG_INFO("returning {}", std::get<1>(*it));
            max = std::get<0>(*it);
            considered_chas_flag[std::get<1>(*it)] = true;
            considered_chas.push_back(std::get<1>(*it));
            break;
        }
        it++;
    }

    // SPDLOG_INFO("communication between threads {} and {} uses the following chas the most", std::get<1>(*it), std::get<0>(*it), max);
    while (it != ranked_cha_access_count_per_pair.end())
    {
        // std::pair<int, int> tid_pair(std::get<2>(*it), std::get<3>(*it));
        if ((std::get<2>(*it) == tid1 && std::get<3>(*it) == tid2) || (std::get<3>(*it) == tid1 && std::get<2>(*it) == tid2))                {
            // SPDLOG_INFO("returning {}", std::get<1>(*it));
            if (considered_chas_flag[std::get<1>(*it)] == false && std::get<0>(*it) > (0.9 * max))
            {
                // SPDLOG_INFO("cha {}, access count: {}, max: {}", std::get<1>(*it), std::get<0>(*it), max);
                considered_chas_flag[std::get<1>(*it)] = true;
                considered_chas.push_back(std::get<1>(*it));
            }
        }
        it++;
    }

    int x_total = 0;
    int y_total = 0;
    int cha_count = 0;
    // SPDLOG_INFO("communication between threads {} and {} involves the following CHAs:");
    for (auto it1 : considered_chas)
    {
        auto tile = topo.getTile(it1);
        x_total += tile.x;
        y_total += tile.y;
        cha_count++;
        // SPDLOG_INFO("cha {}, x: {}, y: {}", tile.cha, tile.x, tile.y);
    }

    assert(cha_count != 0);
    int x_coord = x_total / cha_count;
    int y_coord = y_total / cha_count;
    auto tile = topo.getTile(x_coord, y_coord);
    // SPDLOG_INFO("the center of gravity is cha {}, x: {}, y: {}", tile.cha, tile.x, tile.y);
    // approximate the algorithm now
    return tile.cha;
}

int main (int argc, string argv[])
{
   long c;

   while ((c = getopt(argc, argv, "h")) != -1) {
     switch(c) {
      case 'h':
	Help();
	exit(-1);
	break;
      default:
	fprintf(stderr, "Only valid option is \"-h\".\n");
	exit(-1);
	break;
     }
   }

   Global = NULL;

   pthread_spin_init(&map_spinlock, PTHREAD_PROCESS_SHARED);
   initparam(defv); // modify initparam to read input from stdin only once
   startrun(); // create another version of this function that reuses loaded data
   initoutput(); // no need for modification, can be repeated
   tab_init(); // no need to be recalled in the next iteration of barnes computation

   // the following 5 initializations need to be repeated before each iteration of barnes computation
   Global->tracktime = 0;
   Global->partitiontime = 0;
   Global->treebuildtime = 0;
   Global->forcecalctime = 0;
   Global->current_id = 0;

   using namespace std;
   using std::chrono::duration;
   using std::chrono::duration_cast;
   using std::chrono::high_resolution_clock;
   //using std::chrono::milliseconds;
   using std::chrono::nanoseconds;

#if 0
   uint64_t total_traffic_diff = 0;
    const int traffic_i = (const int) TrafficType::Invalidate;
    const auto traffic_type = static_cast<TrafficType>(traffic_i);
#endif
   CLOCK(Global->computestart);

   printf("COMPUTESTART  = %12lu\n",Global->computestart);

    std::cout << "base cores: ";
    std::vector<int> base_assigned_cores;
    for (int i = 0; i < getCoreCount(); ++i)
    {
        if (i % 2 == 0)
        {
            base_assigned_cores.push_back(i);
            std::cout << i << ' ';
            // this is to bind cores in socket-0. all cores are even numbered in this socket.
        }
    }
    std::cout << std::endl;
    assert(base_assigned_cores.size() == 28); 
//#if 0
   std::cerr << "before SlaveStart\n";
   const auto profiled_base_time_start = high_resolution_clock::now();
   //const auto &before_traffic_vals = storeTraffic(traffic_type);
   CREATE(SlaveStart<true>, static_cast<void*>(base_assigned_cores.data()), NPROC);
//#endif
   WAIT_FOR_END(NPROC);
   const auto profiled_base_time_end = high_resolution_clock::now();
   //const auto &after_traffic_vals = storeTraffic(traffic_type);

   CLOCK(Global->computeend);

   const auto elapsed_profiled_base_time = duration_cast<nanoseconds>(profiled_base_time_end - profiled_base_time_start).count();

   std::cout << "Ended base execution. elapsed time: " << elapsed_profiled_base_time << "ns" << std::endl;
   printf("COMPUTEEND    = %12lu\n",Global->computeend);
   printf("COMPUTETIME   = %12lu\n",Global->computeend - Global->computestart);
   printf("TRACKTIME     = %12lu\n",Global->tracktime);
   printf("PARTITIONTIME = %12lu\t%5.2f\n",Global->partitiontime,
	  ((float)Global->partitiontime)/Global->tracktime);
   printf("TREEBUILDTIME = %12lu\t%5.2f\n",Global->treebuildtime,
	  ((float)Global->treebuildtime)/Global->tracktime);
   printf("FORCECALCTIME = %12lu\t%5.2f\n",Global->forcecalctime,
	  ((float)Global->forcecalctime)/Global->tracktime);
   printf("RESTTIME      = %12lu\t%5.2f\n",
	  Global->tracktime - Global->partitiontime -
	  Global->treebuildtime - Global->forcecalctime,
	  ((float)(Global->tracktime-Global->partitiontime-
		   Global->treebuildtime-Global->forcecalctime))/
	  Global->tracktime);
//#endif
//#if 0   
   // preprocessing begins here
   //std::cout << "Starting preprocesing algo..." << std::endl;
   //const auto algo_start = high_resolution_clock::now();

   assert(NPROC > 1);  // below algo depends on this. we will find thread pairs.

   std::cout << "Starting preprocesing algo..." << std::endl;
   const auto algo_start = high_resolution_clock::now();
#if 0
   auto head = threadid_addresses_map.begin();
   std::cout << "here 1\n";
   std::map<long, std::multiset<double *>>::iterator tail;
   if(head != threadid_addresses_map.end())
   	tail = std::next(threadid_addresses_map.begin());
   std::cout << "here 2\n";

   std::multiset<tuple<int, int, int>, greater<>> total_comm_count_t1_t2;
   std::multiset<tuple<int, int, int, int>, greater<>> total_cha_freq_count_t1_t2;

   // map<pair<int, int>, multiset<Cell *>> pairing_addresses;
   std::cout << "before head != threadid_addresses_map.end()" << std::endl;
   while (head != threadid_addresses_map.end()) {
        const auto orig_tail = tail;
        while (tail != threadid_addresses_map.end()) {
            const int t1 = head->first;
            const int t2 = tail->first;   
	    // cout << "head: " << t1 << ", tail: " << t2 << endl;

            const multiset<double *> t1_addresses = head->second;
            const multiset<double *> t2_addresses = tail->second;

	    std::multiset<double *> common_addresses;
            std::set_intersection(t1_addresses.begin(), t1_addresses.end(), t2_addresses.begin(), t2_addresses.end(),
                                  std::inserter(common_addresses, common_addresses.begin()));

	    std::unordered_map<int, int> cha_freq_map;
	    
	    // this part is changed wrt fluidanimate.
            for(const double* common_addr : common_addresses) {
              ++cha_freq_map[findCha(common_addr)];
            }
            // this part is changed wrt fluidanimate.

	    for(const auto& [cha, freq] : cha_freq_map) {
                total_cha_freq_count_t1_t2.insert({freq, cha, t1, t2});
            }

	    total_comm_count_t1_t2.insert({common_addresses.size(), t1, t2});
	    // pairing_addresses[{t1, t2}] = common_addresses; // pairing is not used at the moment. here just for clarity.
	    
	    ++tail;
	}
	tail = std::next(orig_tail);
	++head;
   }    

   int mapped_thread_count = 0;
   auto it = total_cha_freq_count_t1_t2.begin();
   auto it1 = total_comm_count_t1_t2.begin();
#endif

    //const auto algo_end = high_resolution_clock::now();
    //std::cout << "Ended preprocesing algo. elapsed time: " << duration_cast<milliseconds>(algo_end - algo_start).count() << "ms" << std::endl;
    //std::cout << "Ended preprocesing algo. elapsed time: " << duration_cast<nanoseconds>(algo_end - algo_start).count() << "ns" << std::endl;


    std::multiset<std::tuple<int, int, int, int>, greater<std::tuple<int, int, int, int>>> ranked_cha_access_count_per_pair;
    std::multiset<std::tuple<int, int, int>, greater<std::tuple<int, int, int>>> ranked_communication_count_per_pair;

    for (int i = 0; i < comm_matrix.size(); ++i)
        {
                for (int j = i + 1; j < comm_matrix[i].size(); ++j)
                {
                        if(comm_matrix[i][j].first > 0 || comm_matrix[j][i].first > 0) {
				
                                ranked_communication_count_per_pair.insert({comm_matrix[i][j].first + comm_matrix[j][i].first, i, j});
				for (int k = 0; k < NPROC; k++)
                        	{
					if(comm_matrix[i][j].second[k] > 0 || comm_matrix[j][i].second[k] > 0)
						ranked_cha_access_count_per_pair.insert({comm_matrix[i][j].second[k] + comm_matrix[j][i].second[k], k, i, j});
				}
			}
                }
        }

    int mapped_thread_count = 0;
    auto it = ranked_cha_access_count_per_pair.begin();
    auto it1 = ranked_communication_count_per_pair.begin();
    std::vector<int> thread_to_core(NPROC, -1);
    auto topo = Topology(cha_core_map, CAPID6);
    std::vector<Tile> mapped_tiles;

    while (mapped_tiles.size() < NPROC && it1 != ranked_communication_count_per_pair.end())
        {
                // std::pair<int, int> tid_pair(std::get<2>(*it), std::get<3>(*it));
                std::pair<int, int> tid_pair(std::get<1>(*it1), std::get<2>(*it1));
                if (thread_to_core[tid_pair.first] == -1 || thread_to_core[tid_pair.second] == -1)
                {
                        //SPDLOG_INFO("thread {} and thread {} has {} communication count", tid_pair.first, tid_pair.second, std::get<0>(*it1));
                        int cha_id = getMostAccessedCHA(tid_pair.first, tid_pair.second, ranked_cha_access_count_per_pair, topo);
                        if (cha_id == -1)
                        {
                                //SPDLOG_INFO("error: cha is -1");
                                it1++;
                                continue;
                        }
                        // auto tile = topo.getTile(std::get<1>(*it));
                        auto tile = topo.getTile(cha_id);
                        //SPDLOG_TRACE("cha {}, is colocated with core {}", cha_id, tile.core);
                        if (thread_to_core[tid_pair.first] == -1)
                        {
                                // SPDLOG_INFO("fetching a tile closest to tile with cha {} and core {}, cha supposed to be {}", tile.cha, tile.core, std::get<1>(*it));
                                auto closest_tile = topo.getClosestTile(tile, mapped_tiles);
                                // auto closest_tile = topo.getClosestTilewithThreshold(tile, mapped_tiles);
                                //SPDLOG_TRACE("* closest _available_ core to cha {} is: {}", tile.cha, closest_tile.core);
                                mapped_tiles.push_back(closest_tile);
                                thread_to_core[tid_pair.first] = closest_tile.core;
                                //SPDLOG_TRACE("assigned thread with id {} to core {}", tid_pair.first, closest_tile.core);
                        }

                        if (thread_to_core[tid_pair.second] == -1)
                        {
                                // SPDLOG_INFO("fetching a tile closest to tile with cha {} and core {}, cha supposed to be {}", tile.cha, tile.core, std::get<1>(*it));
                                auto closest_tile = topo.getClosestTile(tile, mapped_tiles);
                                // auto closest_tile = topo.getClosestTilewithThreshold(tile, mapped_tiles);
                                //SPDLOG_TRACE("# closest _available_ core to cha {} is: {}", tile.cha, closest_tile.core);
                                mapped_tiles.push_back(closest_tile);
                                thread_to_core[tid_pair.second] = closest_tile.core;
                                //SPDLOG_TRACE("assigned thread with id {} to core {}", tid_pair.second, closest_tile.core);
                        }
                }

                it1++;
        }
       
	for (int i = 0; i < NPROC; i++) 
        {
                //thread_to_core[i] = i;
                std::vector<int>::iterator it =
                        std::find(thread_to_core.begin(), thread_to_core.end(), i);
                if (it == thread_to_core.end())
                {
                        std::vector<int>::iterator it1 =
                                std::find(thread_to_core.begin(), thread_to_core.end(), -1);
                        *it1 = i;
                }
        }
    const auto algo_end = high_resolution_clock::now();
    //std::cout << "Ended preprocesing algo. elapsed time: " << duration_cast<milliseconds>(algo_end - algo_start).count() << "ms" << std::endl;
    std::cout << "Ended preprocesing algo. elapsed time: " << duration_cast<nanoseconds>(algo_end - algo_start).count() << "ns" << std::endl;

    int ii = 0;
    for (auto ptr : thread_to_core) {
        std::cout << "thread " << ii << " is mapped to core " << ptr << std::endl;
        // SPDLOG_INFO("thread {} is mapped to core {} ", ii, ptr);
        ii++;
    }

    //assert(thread_to_core.size() == NPROC);
    topo.printTopology(); 
//#endif
#if 0

    // base BM.

    //initparam(defv); // modify initparam to read input from stdin only once
    startrun_repeated(); // create another version of this function that reuses loaded data
    initoutput(); // no need for modification, can be repeated
    //tab_init(); // no need to be recalled in the next iteration of barnes computation

    // the following 5 initializations need to be repeated before each iteration of barnes computation
    Global->tracktime = 0;
    Global->partitiontime = 0;
    Global->treebuildtime = 0;
    Global->forcecalctime = 0;
    Global->current_id = 0;

    std::cout << "Now running base BM" << std::endl;
    const auto base_start0 = high_resolution_clock::now();

    CREATE(SlaveStart<false>, static_cast<void*>(base_assigned_cores.data()), NPROC);

    WAIT_FOR_END(NPROC);

    // std::cout << "AFTER JOIN. ended base bm" << std::endl;

    const auto base_end0 = high_resolution_clock::now();
    //const auto elapsed_base = duration_cast<milliseconds>(base_end - base_start).count();
    const auto elapsed_base0 = duration_cast<nanoseconds>(base_end0 - base_start0).count();
    std::cout << "Ended base BM for cache warming. elapsed time: " << elapsed_base0 << "ns" << std::endl; 
#endif
//#endif
#if 0
     //barrier(Global->Barstart,NPROC);

    //initparam(defv); // modify initparam to read input from stdin only once
    startrun_repeated(); // create another version of this function that reuses loaded data
    initoutput(); // no need for modification, can be repeated
    //tab_init(); // no need to be recalled in the next iteration of barnes computation

    // the following 5 initializations need to be repeated before each iteration of barnes computation
    Global->tracktime = 0;
    Global->partitiontime = 0;
    Global->treebuildtime = 0;
    Global->forcecalctime = 0;
    Global->current_id = 0;

    uint64_t total_traffic_diff = 0;
    const int traffic_i = (const int) TrafficType::Invalidate;
    const auto traffic_type = static_cast<TrafficType>(traffic_i);
    // cha aware BM.
    std::cout << "Now running cha aware BM" << std::endl;
    //assert(__threads__<__MAX_THREADS__);
    const auto cha_aware_start = high_resolution_clock::now();

    //const auto &before_traffic_vals = storeTraffic(traffic_type);
    CREATE(SlaveStart<false>, static_cast<void*>(thread_to_core.data() /*base_assigned_cores.data()*/), NPROC);

    WAIT_FOR_END(NPROC); 
    //const auto &after_traffic_vals = storeTraffic(traffic_type);
    // std::cout << "AFTER JOIN. ended cha aware bm" << std::endl;

    const auto cha_aware_end = high_resolution_clock::now();
    //const auto elapsed_cha_aware = duration_cast<milliseconds>(cha_aware_end - cha_aware_start).count();
    const auto elapsed_cha_aware = duration_cast<nanoseconds>(cha_aware_end - cha_aware_start).count();
    std::cout << "Ended cha aware BM. elapsed time: " << elapsed_cha_aware << "ns" << std::endl;
#endif

#if 0
    assert(before_traffic_vals.size() == after_traffic_vals.size());
    for (int i = 0; i < before_traffic_vals.size(); ++i)
    {
	for (int j = 0; j < before_traffic_vals[i].size(); ++j)
	{
		assert(after_traffic_vals[i][j] >= before_traffic_vals[i][j]);
		total_traffic_diff += (after_traffic_vals[i][j] - before_traffic_vals[i][j]);
	}
    } 
    std::cout << "captured traffic: " << total_traffic_diff << std::endl;
#endif
//#if 0
    // base BM.

    //initparam(defv); // modify initparam to read input from stdin only once
    startrun_repeated(); // create another version of this function that reuses loaded data
    initoutput(); // no need for modification, can be repeated
    //tab_init(); // no need to be recalled in the next iteration of barnes computation

    // the following 5 initializations need to be repeated before each iteration of barnes computation
    Global->tracktime = 0;
    Global->partitiontime = 0;
    Global->treebuildtime = 0;
    Global->forcecalctime = 0;
    Global->current_id = 0;

    std::cout << "Now running base BM" << std::endl;
    const auto base_start = high_resolution_clock::now(); 

    CREATE(SlaveStart<false>, static_cast<void*>(base_assigned_cores.data()), NPROC);

    WAIT_FOR_END(NPROC);

    //std::cout << "AFTER JOIN. ended base bm" << std::endl;

    const auto base_end = high_resolution_clock::now();
    //const auto elapsed_base = duration_cast<milliseconds>(base_end - base_start).count();
    const auto elapsed_base = duration_cast<nanoseconds>(base_end - base_start).count();
    std::cout << "Ended base BM. elapsed time: " << elapsed_base << "ns" << std::endl;
  
    //std::cout << "latency improv percentage: " << ((elapsed_base - elapsed_cha_aware) / static_cast<double>(elapsed_base)) * 100 << std::endl;
//#endif
    pthread_spin_destroy(&map_spinlock);
    //std::cerr << "after pthread_spin_destroy\n";
    MAIN_END;
}

/*
 * ANLINIT : initialize ANL macros
 */
void ANLinit()
{
   MAIN_INITENV(,70000000,);
   /* Allocate global, shared memory */

   Global = (struct GlobalMemory *) G_MALLOC(sizeof(struct GlobalMemory));
   if (Global==NULL) error("No initialization for Global\n");

#if 0
   const int ret = posix_memalign((void **)(&Global), CACHELINE_SIZE, sizeof(struct GlobalMemory));
   assert(ret == 0);

   if (Global==NULL) error("No initialization for Global\n"); 
#endif

   // create a variant of barinit that only zeroes out counter and cycle
   BARINIT(Global->Barrier, NPROC); 

   // no need to be added in the variant of ANLinit
   LOCKINIT(Global->CountLock); 

   // no need to be added in the variant of ANLinit
   LOCKINIT(Global->io_lock);
}

/*
 * ANLINIT : modified
 */
void ANLinit_modified()
{
   MAIN_INITENV(,70000000,);
   /* Allocate global, shared memory */

   //Global = (struct GlobalMemory *) G_MALLOC(sizeof(struct GlobalMemory));
   //if (Global==NULL) error("No initialization for Global\n");

   //const int ret = posix_memalign((void **)(&Global), CACHELINE_SIZE, sizeof(struct GlobalMemory));
   //assert(ret == 0);

   //if (Global==NULL) error("No initialization for Global\n");

   // create a variant of barinit that only zeroes out counter and cycle
   BARINIT_NO_INIT(Global->Barrier, NPROC);

   // no need to be added in the variant of ANLinit
   //LOCKINIT(Global->CountLock);

   // no need to be added in the variant of ANLinit
   //LOCKINIT(Global->io_lock);
}

/*
 * INIT_ROOT: Processor 0 reinitialize the global root at each time step
 */
void init_root()
{
   long i;

   Global->G_root=Local[0].ctab;
   Global->G_root->seqnum = 0;
   Type(Global->G_root) = CELL;
   Done(Global->G_root) = FALSE;
   Level(Global->G_root) = IMAX >> 1;
   for (i = 0; i < NSUB; i++) {
      Subp(Global->G_root)[i] = NULL;
   }
   Local[0].mynumcell=1;
}

long Log_base_2(long number)
{
   long cumulative;
   long out;

   cumulative = 1;
   for (out = 0; out < 20; out++) {
      if (cumulative == number) {
         return(out);
      }
      else {
         cumulative = cumulative * 2;
      }
   }

   fprintf(stderr,"Log_base_2: couldn't find log2 of %ld\n", number);
   exit(-1);
}

/*
 * TAB_INIT : allocate body and cell data space
 */

void tab_init()
{
   long i;

   /*allocate leaf/cell space */
   maxleaf = (long) ((double) fleaves * nbody);
   maxcell = fcells * maxleaf;
   for (i = 0; i < NPROC; ++i) {
      //Local[i].ctab = (cellptr) G_MALLOC((maxcell / NPROC) * sizeof(cell));

      const int ret = posix_memalign((void **)(&(Local[i].ctab)), CACHELINE_SIZE, (maxcell / NPROC) * sizeof(cell));
      assert(ret == 0);

      if (Local[i].ctab==NULL) error("No initialization for ctab\n");

      Local[i].ltab = (leafptr) G_MALLOC((maxleaf / NPROC) * sizeof(leaf));
   }

   /*allocate space for personal lists of body pointers */
   maxmybody = (nbody+maxleaf*MAX_BODIES_PER_LEAF)/NPROC;
   Local[0].mybodytab = (bodyptr*) G_MALLOC(NPROC*maxmybody*sizeof(bodyptr));
   /* space is allocated so that every */
   /* process can have a maximum of maxmybody pointers to bodies */
   /* then there is an array of bodies called bodytab which is  */
   /* allocated in the distribution generation or when the distr. */
   /* file is read */
   maxmycell = maxcell / NPROC;
   maxmyleaf = maxleaf / NPROC;
   Local[0].mycelltab = (cellptr*) G_MALLOC(NPROC*maxmycell*sizeof(cellptr));
   Local[0].myleaftab = (leafptr*) G_MALLOC(NPROC*maxmyleaf*sizeof(leafptr));

   //CellLock = (struct CellLockType *) G_MALLOC(sizeof(struct CellLockType));

   const int ret = posix_memalign((void **)(&CellLock), CACHELINE_SIZE, sizeof(struct CellLockType));
   assert(ret == 0);

   if (CellLock==NULL) error("No initialization for CellLock\n");

   ALOCKINIT(CellLock->CL,MAXLOCK);
}

void stick_this_thread_to_core(int core_id) {
    int num_cores = sysconf(_SC_NPROCESSORS_ONLN);
    if (core_id < 0 || core_id >= num_cores) {
        std::cerr << "error binding thread to core: " << core_id << '\n';
        // SPDLOG_ERROR("error binding thread to core {}!", core_id);
        return;
    }

    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(core_id, &cpuset);

    pthread_t current_thread = pthread_self();

    int res = pthread_setaffinity_np(current_thread, sizeof(cpu_set_t), &cpuset);

    if (res == 0) {
        // std::cout << "thread bound to core " << core_id << std::endl;
        //        SPDLOG_INFO("Thread bound to core {} successfully.", core_id);
    } else {
        //        SPDLOG_ERROR("Error in binding this thread to core {}.", core_id);
    }
}

/*
 * SLAVESTART: main task for each processor
 */
template<bool is_preprocessing>
void SlaveStart(void* data)
{
   //printf("SlaveStart begins\n");
   long ProcessId;
   assert(data);
   //std::cerr << "in SlaveStart 1\n";
   int* cores = static_cast<int*>(data);

   /* Get unique ProcessId */
   LOCK(Global->CountLock);
     ProcessId = Global->current_id++;
   UNLOCK(Global->CountLock);

   stick_this_thread_to_core(cores[static_cast<int>(ProcessId)]);

   //std::cerr << "in SlaveStart 2\n";
   BARINCLUDE(Global->Barrier);

/* POSSIBLE ENHANCEMENT:  Here is where one might pin processes to
   processors to avoid migration */

   /* initialize mybodytabs */
   Local[ProcessId].mybodytab = Local[0].mybodytab + (maxmybody * ProcessId);
   /* note that every process has its own copy   */
   /* of mybodytab, which was initialized to the */
   /* beginning of the whole array by proc. 0    */
   /* before create                              */
   Local[ProcessId].mycelltab = Local[0].mycelltab + (maxmycell * ProcessId);
   Local[ProcessId].myleaftab = Local[0].myleaftab + (maxmyleaf * ProcessId);
/* POSSIBLE ENHANCEMENT:  Here is where one might distribute the
   data across physically distributed memories as desired.

   One way to do this is as follows:

   long i;

   if (ProcessId == 0) {
     for (i=0;i<NPROC;i++) {
       Place all addresses x such that
         &(Local[i]) <= x < &(Local[i])+
           sizeof(struct local_memory) on node i
       Place all addresses x such that
         &(Local[i].mybodytab) <= x < &(Local[i].mybodytab)+
           maxmybody * sizeof(bodyptr) - 1 on node i
       Place all addresses x such that
         &(Local[i].mycelltab) <= x < &(Local[i].mycelltab)+
           maxmycell * sizeof(cellptr) - 1 on node i
       Place all addresses x such that
         &(Local[i].myleaftab) <= x < &(Local[i].myleaftab)+
           maxmyleaf * sizeof(leafptr) - 1 on node i
     }
   }

   barrier(Global->Barstart,NPROC);

*/

   Local[ProcessId].tout = Local[0].tout;
   Local[ProcessId].tnow = Local[0].tnow;
   Local[ProcessId].nstep = Local[0].nstep;

   find_my_initial_bodies(bodytab, nbody, ProcessId);

   /* main loop */
   while (Local[ProcessId].tnow < tstop + 0.1 * dtime) {
      stepsystem<is_preprocessing>(ProcessId);
//      printtree(Global->G_root);
      //printf("Going to next step!!!\n");
   }
   //std::cerr << "in SlaveStart 3\n";
   //printf("SlaveStart ends\n");
}


/*
 * STARTRUN: startup hierarchical N-body code.
 */

void startrun()
{
   long seed;
   infile = getparam("in"); // alter getparam to read from the data structure containing data from stdin
   if (*infile != '\0'/*NULL*/) {
      inputdata();
   }
   else {
      nbody = getiparam("nbody");
      std::cerr << "nbody: " << nbody << "\n";
      if (nbody < 1) {
	 error("startrun: absurd nbody\n");
      }
      seed = getiparam("seed");
   }
   outfile = getparam("out");
   dtime = getdparam("dtime");
   dthf = 0.5 * dtime;
   eps = getdparam("eps");
   epssq = eps*eps;
   tol = getdparam("tol");
   tolsq = tol*tol;
   fcells = getdparam("fcells");
   fleaves = getdparam("fleaves");
   tstop = getdparam("tstop");
   dtout = getdparam("dtout");
   NPROC = getiparam("NPROC");
   Local[0].nstep = 0;
   pranset(seed); // no need for modification, can be repeated
   testdata(); // create another variant of this function that does not initialize bodytab and repeat it
   ANLinit(); // create another variant of ANLinit that only calls the variant of barinit
   setbound(); // no need for modification, can be repeated
   Local[0].tout = Local[0].tnow + dtout;
}

void startrun_repeated()
{
   long seed;
   infile = getparam("in"); // alter getparam to read from the data structure containing data from stdin
   if (*infile != '\0'/*NULL*/) {
      inputdata();
   }
   else {
      nbody = getiparam("nbody");
      if (nbody < 1) {
	 error("startrun: absurd nbody\n");
      }
      seed = getiparam("seed");
   }
   outfile = getparam("out");
   dtime = getdparam("dtime");
   dthf = 0.5 * dtime;
   eps = getdparam("eps");
   epssq = eps*eps;
   tol = getdparam("tol");
   tolsq = tol*tol;
   fcells = getdparam("fcells");
   fleaves = getdparam("fleaves");
   tstop = getdparam("tstop");
   dtout = getdparam("dtout");
   NPROC = getiparam("NPROC");
   Local[0].nstep = 0;
   pranset(seed); // no need for modification, can be repeated
   testdata_no_alloc(); // create another variant of this function that does not initialize bodytab and repeat it
   ANLinit_modified(); // create another variant of ANLinit that only calls the variant of barinit
   setbound(); // no need for modification, can be repeated
   Local[0].tout = Local[0].tnow + dtout;
}

/*
 * TESTDATA: generate Plummer model initial conditions for test runs,
 * scaled to units such that M = -4E = G = 1 (Henon, Hegge, etc).
 * See Aarseth, SJ, Henon, M, & Wielen, R (1974) Astr & Ap, 37, 183.
 */

#define MFRAC  0.999                /* mass cut off at MFRAC of total */

void testdata()
{
   real rsc, vsc, r, v, x, y;
   vector cmr, cmv;
   register bodyptr p;
   long rejects = 0;
   long halfnbody, i;
   float offset;
   register bodyptr cp;

   headline = "Hack code: Plummer model";
   Local[0].tnow = 0.0;
   bodytab = (bodyptr) G_MALLOC(nbody * sizeof(body));
   if (bodytab == NULL) {
      error("testdata: not enough memory\n");
   }
   rsc = 9 * PI / 16;
   vsc = sqrt(1.0 / rsc);

   CLRV(cmr);
   CLRV(cmv);

   halfnbody = nbody / 2;
   if (nbody % 2 != 0) halfnbody++;
   for (p = bodytab; p < bodytab+halfnbody; p++) {
      Type(p) = BODY;
      Mass(p) = 1.0 / nbody;
      Cost(p) = 1;

      r = 1 / sqrt(pow(xrand(0.0, MFRAC), -2.0/3.0) - 1);
      /*   reject radii greater than 10 */
      while (r > 9.0) {
	 rejects++;
	 r = 1 / sqrt(pow(xrand(0.0, MFRAC), -2.0/3.0) - 1);
      }
      pickshell(Pos(p), rsc * r);
      ADDV(cmr, cmr, Pos(p));
      do {
	 x = xrand(0.0, 1.0);
	 y = xrand(0.0, 0.1);

      } while (y > x*x * pow(1 - x*x, 3.5));

      v = sqrt(2.0) * x / pow(1 + r*r, 0.25);
      pickshell(Vel(p), vsc * v);
      ADDV(cmv, cmv, Vel(p));
   }

   offset = 4.0;

   for (p = bodytab + halfnbody; p < bodytab+nbody; p++) {
      Type(p) = BODY;
      Mass(p) = 1.0 / nbody;
      Cost(p) = 1;

      cp = p - halfnbody;
      for (i = 0; i < NDIM; i++){
	 Pos(p)[i] = Pos(cp)[i] + offset;
	 Vel(p)[i] = Vel(cp)[i];
      }
      ADDV(cmr, cmr, Pos(p));
      ADDV(cmv, cmv, Vel(p));
   }

   DIVVS(cmr, cmr, (real) nbody);
   DIVVS(cmv, cmv, (real) nbody);

   for (p = bodytab; p < bodytab+nbody; p++) {
      SUBV(Pos(p), Pos(p), cmr);
      SUBV(Vel(p), Vel(p), cmv);
   }
}

void testdata_no_alloc()
{
   real rsc, vsc, r, v, x, y;
   vector cmr, cmv;
   register bodyptr p;
   long rejects = 0;
   long halfnbody, i;
   float offset;
   register bodyptr cp;

   headline = "Hack code: Plummer model";
   Local[0].tnow = 0.0;
   //bodytab = (bodyptr) G_MALLOC(nbody * sizeof(body));
   if (bodytab == NULL) {
      error("testdata: not enough memory\n");
   }
   rsc = 9 * PI / 16;
   vsc = sqrt(1.0 / rsc);

   CLRV(cmr);
   CLRV(cmv);

   halfnbody = nbody / 2;
   if (nbody % 2 != 0) halfnbody++;
   for (p = bodytab; p < bodytab+halfnbody; p++) {
      Type(p) = BODY;
      Mass(p) = 1.0 / nbody;
      Cost(p) = 1;

      r = 1 / sqrt(pow(xrand(0.0, MFRAC), -2.0/3.0) - 1);
      /*   reject radii greater than 10 */
      while (r > 9.0) {
	 rejects++;
	 r = 1 / sqrt(pow(xrand(0.0, MFRAC), -2.0/3.0) - 1);
      }
      pickshell(Pos(p), rsc * r);
      ADDV(cmr, cmr, Pos(p));
      do {
	 x = xrand(0.0, 1.0);
	 y = xrand(0.0, 0.1);

      } while (y > x*x * pow(1 - x*x, 3.5));

      v = sqrt(2.0) * x / pow(1 + r*r, 0.25);
      pickshell(Vel(p), vsc * v);
      ADDV(cmv, cmv, Vel(p));
   }

   offset = 4.0;

   for (p = bodytab + halfnbody; p < bodytab+nbody; p++) {
      Type(p) = BODY;
      Mass(p) = 1.0 / nbody;
      Cost(p) = 1;

      cp = p - halfnbody;
      for (i = 0; i < NDIM; i++){
	 Pos(p)[i] = Pos(cp)[i] + offset;
	 Vel(p)[i] = Vel(cp)[i];
      }
      ADDV(cmr, cmr, Pos(p));
      ADDV(cmv, cmv, Vel(p));
   }

   DIVVS(cmr, cmr, (real) nbody);
   DIVVS(cmv, cmv, (real) nbody);

   for (p = bodytab; p < bodytab+nbody; p++) {
      SUBV(Pos(p), Pos(p), cmr);
      SUBV(Vel(p), Vel(p), cmv);
   }
}

/*
 * PICKSHELL: pick a random point on a sphere of specified radius.
 */

void pickshell(real vec[], real rad)
{
   register long k;
   double rsq, rsc;

   do {
      for (k = 0; k < NDIM; k++) {
	 vec[k] = xrand(-1.0, 1.0);
      }
      DOTVP(rsq, vec, vec);
   } while (rsq > 1.0);

   rsc = rad / sqrt(rsq);
   MULVS(vec, vec, rsc);
}



long intpow(long i, long j)
{
    long k;
    long temp = 1;

    for (k = 0; k < j; k++)
        temp = temp*i;
    return temp;
}

/*
 * MAKETREE: initialize tree structure for hack force calculation.
 */
template<bool is_preprocessing>
void maketree(long ProcessId)
{
   bodyptr p, *pp;

   Local[ProcessId].myncell = 0;
   Local[ProcessId].mynleaf = 0;
   if (ProcessId == 0) {
      Local[ProcessId].mycelltab[Local[ProcessId].myncell++] = Global->G_root;
   }
   Local[ProcessId].Current_Root = (nodeptr) Global->G_root;
   for (pp = Local[ProcessId].mybodytab;
        pp < Local[ProcessId].mybodytab+Local[ProcessId].mynbody; pp++) {
      p = *pp;
      if (Mass(p) != 0.0) {
         Local[ProcessId].Current_Root
            = (nodeptr) loadtree<is_preprocessing>(p, (cellptr) Local[ProcessId].Current_Root,
                                 ProcessId);
      }
      else {
         LOCK(Global->io_lock);
         fprintf(stderr, "Process %ld found body %ld to have zero mass\n",
                 ProcessId, (long) p);
         UNLOCK(Global->io_lock);
      }
   }

   {
        unsigned long   Error, Cycle;
        long            Cancel, Temp;

        Error = pthread_mutex_lock(&(Global->Barrier).mutex);
        if (Error != 0) {
                printf("Error while trying to get lock in barrier.\n");
                exit(-1);
        }

        Cycle = (Global->Barrier).cycle;
        if (++(Global->Barrier).counter != (NPROC)) {
                pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, (int *) &Cancel);
                while (Cycle == (Global->Barrier).cycle) {
                        Error = pthread_cond_wait(&(Global->Barrier).cv, &(Global->Barrier).mutex);
                        if (Error != 0) {
                                break;
                        }
                }
                pthread_setcancelstate(Cancel, (int *) &Temp);
        } else {
                (Global->Barrier).cycle = !(Global->Barrier).cycle;
                (Global->Barrier).counter = 0;
                Error = pthread_cond_broadcast(&(Global->Barrier).cv);
        }
        pthread_mutex_unlock(&(Global->Barrier).mutex);
   }

   hackcofm<is_preprocessing>(ProcessId );
   {
        unsigned long   Error, Cycle;
        long            Cancel, Temp;

        Error = pthread_mutex_lock(&(Global->Barrier).mutex);
        if (Error != 0) {
                printf("Error while trying to get lock in barrier.\n");
                exit(-1);
        }

        Cycle = (Global->Barrier).cycle;
        if (++(Global->Barrier).counter != (NPROC)) {
                pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, (int *) &Cancel);
                while (Cycle == (Global->Barrier).cycle) {
                        Error = pthread_cond_wait(&(Global->Barrier).cv, &(Global->Barrier).mutex);
                        if (Error != 0) {
                                break;
                        }
                }
                pthread_setcancelstate(Cancel, (int *) &Temp);
        } else {
                (Global->Barrier).cycle = !(Global->Barrier).cycle;
                (Global->Barrier).counter = 0;
                Error = pthread_cond_broadcast(&(Global->Barrier).cv);
        }
        pthread_mutex_unlock(&(Global->Barrier).mutex);
   }
}

/*
 *  * LOADTREE: descend tree and insert particle.
 *   */
template<bool is_preprocessing>
nodeptr loadtree(bodyptr p, cellptr root, long ProcessId)
{
   long l, xp[NDIM], xor_arr[NDIM], flag;
   long i, j, root_level;
   bool valid_root;
   long kidIndex;
   volatile nodeptr *volatile qptr, mynode;
   leafptr le;

   intcoord(xp, Pos(p));
   valid_root = TRUE;
   for (i = 0; i < NDIM; i++) {
      xor_arr[i] = xp[i] ^ Local[ProcessId].Root_Coords[i];
   }
   for (i = IMAX >> 1; i > Level(root); i >>= 1) {
      for (j = 0; j < NDIM; j++) {
         if (xor_arr[j] & i) {
            valid_root = FALSE;
            break;
         }
      }
      if (!valid_root) {
         break;
      }
   }
   if (!valid_root) {
      if (root != Global->G_root) {
         root_level = Level(root);
         for (j = i; j > root_level; j >>= 1) {
            root = (cellptr) Parent(root);
         }
         valid_root = TRUE;
         for (i = IMAX >> 1; i > Level(root); i >>= 1) {
            for (j = 0; j < NDIM; j++) {
               if (xor_arr[j] & i) {
                  valid_root = FALSE;
                  break;
               }
            }
            if (!valid_root) {
               printf("P%ld body %ld\n", ProcessId, p - bodytab);
               root = Global->G_root;
            }
         }
      }
   }
   root = Global->G_root;
   mynode = (nodeptr) root;
   kidIndex = subindex(xp, Level(mynode));
   qptr = &Subp(mynode)[kidIndex];

   l = Level(mynode) >> 1;
   flag = TRUE;
   while (flag) {                           /* loop descending tree     */
      if (l == 0) {
         error("not enough levels in tree\n");
      }
      if (*qptr == NULL) {
         /* lock the parent cell */
         ALOCK(CellLock->CL, ((cellptr) mynode)->seqnum % MAXLOCK);
	 if constexpr (is_preprocessing)
      	 { 
        	pthread_spin_lock(&map_spinlock);
#if 0
        	threadid_addresses_map[ProcessId].insert(
          		reinterpret_cast<double*>(
            			reinterpret_cast<uintptr_t>(&(CellLock->CL)) &
            				~(CACHELINE_SIZE - 1)
          		)
        	);
#endif
		inc_comm(ProcessId, reinterpret_cast<uintptr_t>(&(CellLock->CL)) & ~(CACHELINE_SIZE - 1));

#if 0
                threadid_addresses_map[ProcessId].insert(
                        reinterpret_cast<double*>(
                                reinterpret_cast<uintptr_t>(&(((cellptr) mynode)->seqnum)) &
                                        ~(CACHELINE_SIZE - 1)
                        )
                );
#endif
		inc_comm(ProcessId, reinterpret_cast<uintptr_t>(&(((cellptr) mynode)->seqnum)) & ~(CACHELINE_SIZE - 1));
                pthread_spin_unlock(&map_spinlock);
      	}
         if (*qptr == NULL) {
            le = InitLeaf((cellptr) mynode, ProcessId);
            Parent(p) = (nodeptr) le;
            Level(p) = l;
            ChildNum(p) = le->num_bodies;
            ChildNum(le) = kidIndex;
            Bodyp(le)[le->num_bodies++] = p;
            *qptr = (nodeptr) le;
            flag = FALSE;
         }
         AULOCK(CellLock->CL, ((cellptr) mynode)->seqnum % MAXLOCK);
	 if constexpr (is_preprocessing)
         {
                pthread_spin_lock(&map_spinlock);

#if 0
                threadid_addresses_map[ProcessId].insert(
                        reinterpret_cast<double*>(
                                reinterpret_cast<uintptr_t>(&(CellLock->CL)) &
                                        ~(CACHELINE_SIZE - 1)
                        )
                );
#endif
		inc_comm(ProcessId, reinterpret_cast<uintptr_t>(&(CellLock->CL)) & ~(CACHELINE_SIZE - 1));
#if 0
		threadid_addresses_map[ProcessId].insert(
                        reinterpret_cast<double*>(
                                reinterpret_cast<uintptr_t>(&(((cellptr) mynode)->seqnum)) &
                                        ~(CACHELINE_SIZE - 1)
                        )
                );
#endif
		inc_comm(ProcessId, reinterpret_cast<uintptr_t>(&(((cellptr) mynode)->seqnum)) & ~(CACHELINE_SIZE - 1));		
                pthread_spin_unlock(&map_spinlock);
        }
         /* unlock the parent cell */
      }
      if (flag && *qptr && (Type(*qptr) == LEAF)) {
         /*   reached a "leaf"?      */
         ALOCK(CellLock->CL, ((cellptr) mynode)->seqnum % MAXLOCK);
	 if constexpr (is_preprocessing)
         {
                pthread_spin_lock(&map_spinlock);

#if 0
                threadid_addresses_map[ProcessId].insert(
                        reinterpret_cast<double*>(
                                reinterpret_cast<uintptr_t>(&(CellLock->CL)) &
                                        ~(CACHELINE_SIZE - 1)
                        )
                );
#endif
		inc_comm(ProcessId, reinterpret_cast<uintptr_t>(&(CellLock->CL)) & ~(CACHELINE_SIZE - 1));
#if 0
		threadid_addresses_map[ProcessId].insert(
                        reinterpret_cast<double*>(
                                reinterpret_cast<uintptr_t>(&(((cellptr) mynode)->seqnum)) &
                                        ~(CACHELINE_SIZE - 1)
                        )
                );
#endif
		inc_comm(ProcessId, reinterpret_cast<uintptr_t>(&(((cellptr) mynode)->seqnum)) & ~(CACHELINE_SIZE - 1));
                pthread_spin_unlock(&map_spinlock);
         }
         /* lock the parent cell */
         if (Type(*qptr) == LEAF) {             /* still a "leaf"?      */
            le = (leafptr) *qptr;
            if (le->num_bodies == MAX_BODIES_PER_LEAF) {
               *qptr = (nodeptr) SubdivideLeaf<is_preprocessing>(le, (cellptr) mynode, l,
                                                  ProcessId);
            }
            else {
               Parent(p) = (nodeptr) le;
               Level(p) = l;
               ChildNum(p) = le->num_bodies;
               Bodyp(le)[le->num_bodies++] = p;
               flag = FALSE;
            }
         }
         AULOCK(CellLock->CL, ((cellptr) mynode)->seqnum % MAXLOCK);
	 if constexpr (is_preprocessing)
         {
                pthread_spin_lock(&map_spinlock);

#if 0
                threadid_addresses_map[ProcessId].insert(
                        reinterpret_cast<double*>(
                                reinterpret_cast<uintptr_t>(&(CellLock->CL)) &
                                        ~(CACHELINE_SIZE - 1)
                        )
                );
#endif
		inc_comm(ProcessId, reinterpret_cast<uintptr_t>(&(CellLock->CL)) & ~(CACHELINE_SIZE - 1));
#if 0
		threadid_addresses_map[ProcessId].insert(
                        reinterpret_cast<double*>(
                                reinterpret_cast<uintptr_t>(&(((cellptr) mynode)->seqnum)) &
                                        ~(CACHELINE_SIZE - 1)
                        )
                );
#endif
		inc_comm(ProcessId, reinterpret_cast<uintptr_t>(&(((cellptr) mynode)->seqnum)) & ~(CACHELINE_SIZE - 1));
                pthread_spin_unlock(&map_spinlock);
        }
         /* unlock the node           */
      }
      if (flag) {
         mynode = *qptr;
         kidIndex = subindex(xp, l);
         qptr = &Subp(*qptr)[kidIndex];  /* move down one level  */
         l = l >> 1;                            /* and test next bit    */
      }
   }
   SETV(Local[ProcessId].Root_Coords, xp);
   return Parent((leafptr) *qptr);
}

template<bool is_preprocessing>
cellptr SubdivideLeaf(leafptr le, cellptr parent, long l, long ProcessId)
{  
   cellptr c;
   long i, index;
   long xp[NDIM];
   bodyptr bodies[MAX_BODIES_PER_LEAF];
   long num_bodies;
   bodyptr p;
   
   /* first copy leaf's bodies to temp array, so we can reuse the leaf */
   num_bodies = le->num_bodies;
   for (i = 0; i < num_bodies; i++) {
      bodies[i] = Bodyp(le)[i];
      Bodyp(le)[i] = NULL;
   }
   le->num_bodies = 0;
   /* create the parent cell for this subtree */
   c = InitCell<is_preprocessing>(parent, ProcessId);
   ChildNum(c) = ChildNum(le);
   /* do first particle separately, so we can reuse le */
   p = bodies[0];
   intcoord(xp, Pos(p));
   index = subindex(xp, l);
   Subp(c)[index] = (nodeptr) le;
   ChildNum(le) = index;
   Parent(le) = (nodeptr) c;
   Level(le) = l >> 1;
   /* set stuff for body */
   Parent(p) = (nodeptr) le;
   ChildNum(p) = le->num_bodies;
   Level(p) = l >> 1;
   /* insert the body */
   Bodyp(le)[le->num_bodies++] = p;
   /* now handle the rest */
   for (i = 1; i < num_bodies; i++) {
      p = bodies[i];
      intcoord(xp, Pos(p));
      index = subindex(xp, l);
      if (!Subp(c)[index]) {
         le = InitLeaf(c, ProcessId);
         ChildNum(le) = index;
         Subp(c)[index] = (nodeptr) le;
      }
      else {
         le = (leafptr) Subp(c)[index];
      }
      Parent(p) = (nodeptr) le;
      ChildNum(p) = le->num_bodies;
      Level(p) = l >> 1;
      Bodyp(le)[le->num_bodies++] = p;
   }
   return c;
}


template<bool is_preprocessing>
cellptr InitCell(cellptr parent, long ProcessId)
{
   cellptr c;

   c = makecell<is_preprocessing>(ProcessId);
   c->processor = ProcessId;
   c->next = NULL;
   c->prev = NULL;
   if (parent == NULL)
      Level(c) = IMAX >> 1;
   else
      Level(c) = Level(parent) >> 1;
   Parent(c) = (nodeptr) parent;
   ChildNum(c) = 0;
   return (c);
}


/*
 *  * MAKECELL: allocation routine for cells.
 *   */
template<bool is_preprocessing>
cellptr makecell(long ProcessId)
{
   cellptr c;
   long i, Mycell;

   if (Local[ProcessId].mynumcell == maxmycell) {
      error("makecell: Proc %ld needs more than %ld cells; increase fcells\n",
            ProcessId,maxmycell);
   }
   Mycell = Local[ProcessId].mynumcell++;
   c = Local[ProcessId].ctab + Mycell;
   c->seqnum = ProcessId*maxmycell+Mycell;
   if constexpr (is_preprocessing)
         {
                pthread_spin_lock(&map_spinlock);

#if 0
                threadid_addresses_map[ProcessId].insert(
                        reinterpret_cast<double*>(
                                reinterpret_cast<uintptr_t>(&(c->seqnum)) &
                                        ~(CACHELINE_SIZE - 1)
                        )
                );
#endif
		inc_comm(ProcessId, reinterpret_cast<uintptr_t>(&(c->seqnum)) & ~(CACHELINE_SIZE - 1));
                pthread_spin_unlock(&map_spinlock);
        } 
   Type(c) = CELL;
   Done(c) = FALSE;
   Mass(c) = 0.0;
   for (i = 0; i < NSUB; i++) {
      Subp(c)[i] = NULL;
   }
   Local[ProcessId].mycelltab[Local[ProcessId].myncell++] = c;
   return (c);
}

/*
 *  * HACKCOFM: descend tree finding center-of-mass coordinates.
 *   */
template<bool is_preprocessing>
void hackcofm(long ProcessId)
{
   long i;
   nodeptr r;
   leafptr l;
   leafptr* ll;
   bodyptr p;
   cellptr q;
   cellptr *cc;
   vector tmpv;

   /* get a cell using get*sub.  Cells are got in reverse of the order in */
   /* the cell array; i.e. reverse of the order in which they were created */
   /* this way, we look at child cells before parents                    */

   for (ll = Local[ProcessId].myleaftab + Local[ProcessId].mynleaf - 1;
        ll >= Local[ProcessId].myleaftab; ll--) {
      l = *ll;
      Mass(l) = 0.0;
      Cost(l) = 0;
      CLRV(Pos(l));
      for (i = 0; i < l->num_bodies; i++) {
         p = Bodyp(l)[i];
         Mass(l) += Mass(p);
         Cost(l) += Cost(p);
         MULVS(tmpv, Pos(p), Mass(p));
         ADDV(Pos(l), Pos(l), tmpv);
      }
      DIVVS(Pos(l), Pos(l), Mass(l));
#ifdef QUADPOLE
      CLRM(Quad(l));
      for (i = 0; i < l->num_bodies; i++) {
         p = Bodyp(l)[i];
         SUBV(dr, Pos(p), Pos(l));
         OUTVP(drdr, dr, dr);
         DOTVP(drsq, dr, dr);
         SETMI(Idrsq);
         MULMS(Idrsq, Idrsq, drsq);
         MULMS(tmpm, drdr, 3.0);
         SUBM(tmpm, tmpm, Idrsq);
         MULMS(tmpm, tmpm, Mass(p));
         ADDM(Quad(l), Quad(l), tmpm);
      }
#endif
      Done(l)=TRUE;
   }
   for (cc = Local[ProcessId].mycelltab+Local[ProcessId].myncell-1;
        cc >= Local[ProcessId].mycelltab; cc--) {
      q = *cc;
      Mass(q) = 0.0;
      Cost(q) = 0;
      CLRV(Pos(q));
      for (i = 0; i < NSUB; i++) {
         r = Subp(q)[i];
         if (r != NULL) {
            while(!Done(r)) {
               /* wait */
            }
            Mass(q) += Mass(r);
            Cost(q) += Cost(r);
            MULVS(tmpv, Pos(r), Mass(r));
            ADDV(Pos(q), Pos(q), tmpv);

	    if constexpr (is_preprocessing)
            {
                pthread_spin_lock(&map_spinlock);
#if 0
                threadid_addresses_map[ProcessId].insert(
                        reinterpret_cast<double*>(
                                reinterpret_cast<uintptr_t>(&(Pos(q))) &
                                        ~(CACHELINE_SIZE - 1)
                        )
                );
#endif
		inc_comm(ProcessId, reinterpret_cast<uintptr_t>(&(Pos(q))) & ~(CACHELINE_SIZE - 1));
                pthread_spin_unlock(&map_spinlock);
           }

            Done(r) = FALSE;
         }
      }
      DIVVS(Pos(q), Pos(q), Mass(q));
#ifdef QUADPOLE
      CLRM(Quad(q));
      for (i = 0; i < NSUB; i++) {
         r = Subp(q)[i];
         if (r != NULL) {
            SUBV(dr, Pos(r), Pos(q));
            OUTVP(drdr, dr, dr);
            DOTVP(drsq, dr, dr);
            SETMI(Idrsq);
            MULMS(Idrsq, Idrsq, drsq);
            MULMS(tmpm, drdr, 3.0);
            SUBM(tmpm, tmpm, Idrsq);
            MULMS(tmpm, tmpm, Mass(r));
            ADDM(tmpm, tmpm, Quad(r));
            ADDM(Quad(q), Quad(q), tmpm);
         }
      }
#endif
      Done(q)=TRUE;
   }
}

/*
 * STEPSYSTEM: advance N-body system one time-step.
 */
template<bool is_preprocessing>
void stepsystem(long ProcessId)
{
    long i;
    real Cavg;
    bodyptr p,*pp;
    vector dvel, vel1, dpos;
    long trackstart, trackend;
    long partitionstart, partitionend;
    long treebuildstart, treebuildend;
    long forcecalcstart, forcecalcend;

    if (Local[ProcessId].nstep == 2) {
/* POSSIBLE ENHANCEMENT:  Here is where one might reset the
   statistics that one is measuring about the parallel execution */
    }

    if ((ProcessId == 0) && (Local[ProcessId].nstep >= 2)) {
        CLOCK(trackstart);
    }

    if (ProcessId == 0) {
       init_root();
    }
    else {
       Local[ProcessId].mynumcell = 0;
       Local[ProcessId].mynumleaf = 0;
    }


    /* start at same time */

    {
    unsigned long   Error, Cycle;
        long            Cancel, Temp;

        Error = pthread_mutex_lock(&(Global->Barrier).mutex);
        if (Error != 0) {
                printf("Error while trying to get lock in barrier.\n");
                exit(-1);
        }

        Cycle = (Global->Barrier).cycle;
        if (++(Global->Barrier).counter != (NPROC)) {
                pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, (int *) &Cancel);
                while (Cycle == (Global->Barrier).cycle) {
                        Error = pthread_cond_wait(&(Global->Barrier).cv, &(Global->Barrier).mutex);
                        if (Error != 0) {
                                break;
                        }
                }
                pthread_setcancelstate(Cancel, (int *) &Temp);
        } else {
                (Global->Barrier).cycle = !(Global->Barrier).cycle;
                (Global->Barrier).counter = 0;
                Error = pthread_cond_broadcast(&(Global->Barrier).cv);
        }
        pthread_mutex_unlock(&(Global->Barrier).mutex);
    }

    if ((ProcessId == 0) && (Local[ProcessId].nstep >= 2)) {
        CLOCK(treebuildstart);
    }

    /* load bodies into tree   */
    maketree<is_preprocessing>(ProcessId);
    if ((ProcessId == 0) && (Local[ProcessId].nstep >= 2)) {
        CLOCK(treebuildend);
        Global->treebuildtime += treebuildend - treebuildstart;
    }

    Housekeep(ProcessId);

    Cavg = (real) Cost(Global->G_root) / (real)NPROC ;
    Local[ProcessId].workMin = (long) (Cavg * ProcessId);
    Local[ProcessId].workMax = (long) (Cavg * (ProcessId + 1)
				      + (ProcessId == (NPROC - 1)));

    if ((ProcessId == 0) && (Local[ProcessId].nstep >= 2)) {
        CLOCK(partitionstart);
    }

    Local[ProcessId].mynbody = 0;
    find_my_bodies(reinterpret_cast<nodeptr>(Global->G_root), 0, BRC_FUC, ProcessId );

/*     B*RRIER(Global->Barcom,NPROC); */
    if ((ProcessId == 0) && (Local[ProcessId].nstep >= 2)) {
        CLOCK(partitionend);
        Global->partitiontime += partitionend - partitionstart;
    }

    if ((ProcessId == 0) && (Local[ProcessId].nstep >= 2)) {
        CLOCK(forcecalcstart);
    }

    ComputeForces<is_preprocessing>(ProcessId);

    if ((ProcessId == 0) && (Local[ProcessId].nstep >= 2)) {
        CLOCK(forcecalcend);
        Global->forcecalctime += forcecalcend - forcecalcstart;
    }

    /* advance my bodies */
    for (pp = Local[ProcessId].mybodytab;
	 pp < Local[ProcessId].mybodytab+Local[ProcessId].mynbody; pp++) {
       p = *pp;
       MULVS(dvel, Acc(p), dthf);
       ADDV(vel1, Vel(p), dvel);
       MULVS(dpos, vel1, dtime);
       ADDV(Pos(p), Pos(p), dpos);
       ADDV(Vel(p), vel1, dvel);

       for (i = 0; i < NDIM; i++) {
          if (Pos(p)[i]<Local[ProcessId].min[i]) {
	     Local[ProcessId].min[i]=Pos(p)[i];
	  }
          if (Pos(p)[i]>Local[ProcessId].max[i]) {
	     Local[ProcessId].max[i]=Pos(p)[i] ;
	  }
       }
    }
    LOCK(Global->CountLock);
    for (i = 0; i < NDIM; i++) {
       if (Global->min[i] > Local[ProcessId].min[i]) {
	  Global->min[i] = Local[ProcessId].min[i];
       }
       if (Global->max[i] < Local[ProcessId].max[i]) {
	  Global->max[i] = Local[ProcessId].max[i];
       }
    }
    UNLOCK(Global->CountLock);

    /* bar needed to make sure that every process has computed its min */
    /* and max coordinates, and has accumulated them into the global   */
    /* min and max, before the new dimensions are computed	       */

    {
        unsigned long   Error, Cycle;
        long            Cancel, Temp;

        Error = pthread_mutex_lock(&(Global->Barrier).mutex);
        if (Error != 0) {
                printf("Error while trying to get lock in barrier.\n");
                exit(-1);
        }

        Cycle = (Global->Barrier).cycle;
        if (++(Global->Barrier).counter != (NPROC)) {
                pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, (int *) &Cancel);
                while (Cycle == (Global->Barrier).cycle) {
                        Error = pthread_cond_wait(&(Global->Barrier).cv, &(Global->Barrier).mutex);
                        if (Error != 0) {
                                break;
                        }
                }
                pthread_setcancelstate(Cancel, (int *) &Temp);
        } else {
                (Global->Barrier).cycle = !(Global->Barrier).cycle;
                (Global->Barrier).counter = 0;
                Error = pthread_cond_broadcast(&(Global->Barrier).cv);
        }
        pthread_mutex_unlock(&(Global->Barrier).mutex);
    }

    if ((ProcessId == 0) && (Local[ProcessId].nstep >= 2)) {
        CLOCK(trackend);
        Global->tracktime += trackend - trackstart;
    }
    if (ProcessId==0) {
      Global->rsize=0;
      SUBV(Global->max,Global->max,Global->min);
      for (i = 0; i < NDIM; i++) {
	if (Global->rsize < Global->max[i]) {
	   Global->rsize = Global->max[i];
	}
      }
      ADDVS(Global->rmin,Global->min,-Global->rsize/100000.0);
      Global->rsize = 1.00002*Global->rsize;
      SETVS(Global->min,1E99);
      SETVS(Global->max,-1E99);
    }
    Local[ProcessId].nstep++;
    Local[ProcessId].tnow = Local[ProcessId].tnow + dtime;
}


/*
 *  * HACKGRAV: evaluate grav field at a given particle.
 *   */
template<bool is_preprocessing>
void hackgrav(bodyptr p, long ProcessId)
{
   Local[ProcessId].pskip = p;
   SETV(Local[ProcessId].pos0, Pos(p));
   Local[ProcessId].phi0 = 0.0;
   CLRV(Local[ProcessId].acc0);
   Local[ProcessId].myn2bterm = 0;
   Local[ProcessId].mynbcterm = 0;
   Local[ProcessId].skipself = FALSE;
   hackwalk<is_preprocessing>(ProcessId);
   Phi(p) = Local[ProcessId].phi0;
   SETV(Acc(p), Local[ProcessId].acc0);
#ifdef QUADPOLE
   Cost(p) = Local[ProcessId].myn2bterm + NDIM * Local[ProcessId].mynbcterm;
#else
   Cost(p) = Local[ProcessId].myn2bterm + Local[ProcessId].mynbcterm;
#endif
}

/*
 *  * HACKWALK: walk the tree opening cells too close to a given point.
 *   */
template<bool is_preprocessing>
void hackwalk(long ProcessId)
{
    walksub<is_preprocessing>(reinterpret_cast <nodeptr>(Global->G_root), Global->rsize * Global->rsize, ProcessId);
}

/*
 *  * WALKSUB: recursive routine to do hackwalk operation.
 *   */
template<bool is_preprocessing>
void walksub(nodeptr n, real dsq, long ProcessId)
{
   nodeptr* nn;
   leafptr l;
   bodyptr p;
   long i;

   if (subdivp<is_preprocessing>(n, dsq, ProcessId)) {
      if (Type(n) == CELL) {
         for (nn = Subp(n); nn < Subp(n) + NSUB; nn++) {
            if (*nn != NULL) {
               walksub<is_preprocessing>(*nn, dsq / 4.0, ProcessId);
            }
         }
      }
      else {
         l = (leafptr) n;
         for (i = 0; i < l->num_bodies; i++) {
            p = Bodyp(l)[i];
            if (p != Local[ProcessId].pskip) {
               gravsub(reinterpret_cast <nodeptr>(p), ProcessId);
            }
            else {
               Local[ProcessId].skipself = TRUE;
            }
         }
      }
   }
   else {
      gravsub(n, ProcessId);
   }
}

/*
 *  * SUBDIVP: decide if a node should be opened.
 *   * Side effects: sets  pmem,dr, and drsq.
 *    */
template<bool is_preprocessing>
bool subdivp(register nodeptr p, real dsq, long ProcessId)
{
   SUBV(Local[ProcessId].dr, Pos(p), Local[ProcessId].pos0);

   if constexpr (is_preprocessing)
   {
   	pthread_spin_lock(&map_spinlock);
#if 0
   	threadid_addresses_map[ProcessId].insert(
        reinterpret_cast<double*>(
        	reinterpret_cast<uintptr_t>(&(Pos(p))) &
             		~(CACHELINE_SIZE - 1)
                )
        );
#endif
	inc_comm(ProcessId, reinterpret_cast<uintptr_t>(&(Pos(p))) & ~(CACHELINE_SIZE - 1));
	pthread_spin_unlock(&map_spinlock);
   }
   DOTVP(Local[ProcessId].drsq, Local[ProcessId].dr, Local[ProcessId].dr);
   Local[ProcessId].pmem = p;
   return (tolsq * Local[ProcessId].drsq < dsq);
}

template<bool is_preprocessing>
void ComputeForces(long ProcessId)
{
   bodyptr p,*pp;
   vector acc1, dacc, dvel;

   for (pp = Local[ProcessId].mybodytab;
	pp < Local[ProcessId].mybodytab+Local[ProcessId].mynbody;pp++) {
      p = *pp;
      SETV(acc1, Acc(p));
      Cost(p)=0;
      hackgrav<is_preprocessing>(p,ProcessId);
      Local[ProcessId].myn2bcalc += Local[ProcessId].myn2bterm;
      Local[ProcessId].mynbccalc += Local[ProcessId].mynbcterm;
      if (!Local[ProcessId].skipself) {       /*   did we miss self-int?  */
	 Local[ProcessId].myselfint++;        /*   count another goofup   */
      }
      if (Local[ProcessId].nstep > 0) {
	 /*   use change in accel to make 2nd order correction to vel      */
	 SUBV(dacc, Acc(p), acc1);
	 MULVS(dvel, dacc, dthf);
	 ADDV(Vel(p), Vel(p), dvel);
      }
   }
}

/*
 * FIND_MY_INITIAL_BODIES: puts into mybodytab the initial list of bodies
 * assigned to the processor.
 */

void find_my_initial_bodies(bodyptr btab, long nbody, long ProcessId)
{
  long extra,offset,i;

  Local[ProcessId].mynbody = nbody / NPROC;
  extra = nbody % NPROC;
  if (ProcessId < extra) {
    Local[ProcessId].mynbody++;
    offset = Local[ProcessId].mynbody * ProcessId;
  }
  if (ProcessId >= extra) {
    offset = (Local[ProcessId].mynbody+1) * extra + (ProcessId - extra)
       * Local[ProcessId].mynbody;
  }
  for (i=0; i < Local[ProcessId].mynbody; i++) {
     Local[ProcessId].mybodytab[i] = &(btab[offset+i]);
  }

  {
        unsigned long   Error, Cycle;
        long            Cancel, Temp;

        Error = pthread_mutex_lock(&(Global->Barrier).mutex);
        if (Error != 0) {
                printf("Error while trying to get lock in barrier.\n");
                exit(-1);
        }

        Cycle = (Global->Barrier).cycle;
        if (++(Global->Barrier).counter != (NPROC)) {
                pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, (int *) &Cancel);
                while (Cycle == (Global->Barrier).cycle) {
                        Error = pthread_cond_wait(&(Global->Barrier).cv, &(Global->Barrier).mutex);
                        if (Error != 0) {
                                break;
                        }
                }
                pthread_setcancelstate(Cancel, (int *) &Temp);
        } else {
                (Global->Barrier).cycle = !(Global->Barrier).cycle;
                (Global->Barrier).counter = 0;
                Error = pthread_cond_broadcast(&(Global->Barrier).cv);
        }
        pthread_mutex_unlock(&(Global->Barrier).mutex);
  } 

}


void find_my_bodies(nodeptr mycell, long work, long direction, long ProcessId)
{
   long i;
   leafptr l;
   nodeptr qptr;

   if (Type(mycell) == LEAF) {
      l = (leafptr) mycell;
      for (i = 0; i < l->num_bodies; i++) {
	 if (work >= Local[ProcessId].workMin - .1) {
	    if((Local[ProcessId].mynbody+2) > maxmybody) {
	       error("find_my_bodies: Processor %ld needs more than %ld bodies; increase fleaves\n", ProcessId, maxmybody);
	    }
	    Local[ProcessId].mybodytab[Local[ProcessId].mynbody++] =
	       Bodyp(l)[i];
	 }
	 work += Cost(Bodyp(l)[i]);
	 if (work >= Local[ProcessId].workMax-.1) {
	    break;
	 }
      }
   }
   else {
      for(i = 0; (i < NSUB) && (work < (Local[ProcessId].workMax - .1)); i++){
	 qptr = Subp(mycell)[Child_Sequence[direction][i]];
	 if (qptr!=NULL) {
	    if ((work+Cost(qptr)) >= (Local[ProcessId].workMin -.1)) {
	       find_my_bodies(qptr,work, Direction_Sequence[direction][i],
			      ProcessId);
	    }
	    work += Cost(qptr);
	 }
      }
   }
}

/*
 * HOUSEKEEP: reinitialize the different variables (in particular global
 * variables) between each time step.
 */

void Housekeep(long ProcessId)
{
   Local[ProcessId].myn2bcalc = Local[ProcessId].mynbccalc
      = Local[ProcessId].myselfint = 0;
   SETVS(Local[ProcessId].min,1E99);
   SETVS(Local[ProcessId].max,-1E99);
}

/*
 * SETBOUND: Compute the initial size of the root of the tree; only done
 * before first time step, and only processor 0 does it
 */
void setbound()
{
   long i;
   real side ;
   bodyptr p;

   SETVS(Local[0].min,1E99);
   SETVS(Local[0].max,-1E99);
   side=0;

   for (p = bodytab; p < bodytab+nbody; p++) {
      for (i=0; i<NDIM;i++) {
	 if (Pos(p)[i]<Local[0].min[i]) Local[0].min[i]=Pos(p)[i] ;
	 if (Pos(p)[i]>Local[0].max[i])  Local[0].max[i]=Pos(p)[i] ;
      }
   }

   SUBV(Local[0].max,Local[0].max,Local[0].min);
   for (i=0; i<NDIM;i++) if (side<Local[0].max[i]) side=Local[0].max[i];
   ADDVS(Global->rmin,Local[0].min,-side/100000.0);
   Global->rsize = 1.00002*side;
   SETVS(Global->max,-1E99);
   SETVS(Global->min,1E99);
}

void Help()
{
   printf("There are a total of twelve parameters, and all of them have default values.\n");
   printf("\n");
   printf("1) infile (char*) : The name of an input file that contains particle data.  \n");
   printf("    The format of the file is:\n");
   printf("\ta) An int representing the number of particles in the distribution\n");
   printf("\tb) An int representing the dimensionality of the problem (3-D)\n");
   printf("\tc) A double representing the current time of the simulation\n");
   printf("\td) Doubles representing the masses of all the particles\n");
   printf("\te) A vector (length equal to the dimensionality) of doubles\n");
   printf("\t   representing the positions of all the particles\n");
   printf("\tf) A vector (length equal to the dimensionality) of doubles\n");
   printf("\t   representing the velocities of all the particles\n");
   printf("\n");
   printf("    Each of these numbers can be separated by any amount of whitespace.\n");
   printf("\n");
   printf("2) nbody (int) : If no input file is specified (the first line is blank), this\n");
   printf("    number specifies the number of particles to generate under a plummer model.\n");
   printf("    Default is 16384.\n");
   printf("\n");
   printf("3) seed (int) : The seed used by the random number generator.\n");
   printf("    Default is 123.\n");
   printf("\n");
   printf("4) outfile (char*) : The name of the file that snapshots will be printed to. \n");
   printf("    This feature has been disabled in the SPLASH release.\n");
   printf("    Default is NULL.\n");
   printf("\n");
   printf("5) dtime (double) : The integration time-step.\n");
   printf("    Default is 0.025.\n");
   printf("\n");
   printf("6) eps (double) : The usual potential softening\n");
   printf("    Default is 0.05.\n");
   printf("\n");
   printf("7) tol (double) : The cell subdivision tolerance.\n");
   printf("    Default is 1.0.\n");
   printf("\n");
   printf("8) fcells (double) : The total number of cells created is equal to \n");
   printf("    fcells * number of leaves.\n");
   printf("    Default is 2.0.\n");
   printf("\n");
   printf("9) fleaves (double) : The total number of leaves created is equal to  \n");
   printf("    fleaves * nbody.\n");
   printf("    Default is 0.5.\n");
   printf("\n");
   printf("10) tstop (double) : The time to stop integration.\n");
   printf("    Default is 0.075.\n");
   printf("\n");
   printf("11) dtout (double) : The data-output interval.\n");
   printf("    Default is 0.25.\n");
   printf("\n");
   printf("12) NPROC (int) : The number of processors.\n");
   printf("    Default is 1.\n");
}
