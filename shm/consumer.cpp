#include "protocol.h"

#include <cstdlib>
#include <sys/mman.h>
#include <cinttypes>
#include <cstring>

#include <iostream>
#include <array>
#include <cassert>
#include <chrono>

using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::milliseconds;
using std::chrono::microseconds;
using std::chrono::nanoseconds;

#if 0
        template <typename T>
std::vector<std::pair<int, int>> findChasOfArrayWithShm(T *arr, int offset, int size)
{
        std::vector<std::pair<int, int>> res(size, {-1, -1});

        for (int i = 0; i < size; ++i) {
    // std::cout << "data[" << i << "]: " << data[i] << ", CHA of &data[" << i << "]: " << data[i + NUM] << std::endl;
                res[i].first = 0; // assumed socket-0.
                res[i].second = (int) arr[offset + i];
        }

        return res;
}
#endif

template <typename T> std::vector<T> findWithShm(const T *const data)
{
    //int array_size = NUM;
    std::vector<T> res(NUM, -1);

  for (int i = 0; i < NUM; ++i) {
    res[i] = (T) data[i + NUM];
    //if((i % 64 == 0) && i < 640)
	//std::cout << "data[" << i << "]: " << (int) data[i] << ", CHA of &data[" << i << "]: " << (int) res[i] << std::endl;
  }

    return res;
}

template <typename T> std::vector<T> findWithHashFunc(const T *const data)
{
    //std::array<T, NUM> res;
    std::vector<T> res(NUM, -1);

  for (int i = 0; i < NUM; ++i) {
    res[i] = findCHAByHashing(uintptr_t(&data[i]), base_sequence_28_skx);
    //if((i % 64 == 0) && i < 640)
        //std::cout << "data[" << i << "]: " << (int) data[i] << ", CHA1 of &data[" << i << "]: " << (int) res[i] << std::endl;
  }

    return res;
}

template <typename T>
std::array<int, NUM> findChasOfArray(T *arr)
{
    std::array<int, NUM> res;

    auto addr = static_cast<T *>(&arr[0]);

    // #pragma omp parallel for
    for (int i = 0; i < NUM;)
    {
        addr = static_cast<T *>(&arr[i]);
        if (addr)
        {
            // SPDLOG_INFO("printing the int the addr points to in order to make it page allocated! => {}", *addr);
            // const void* addr2 = addr; /// important to do this, so that a page would be mapped for the address by OS. this is _important_ for finding cha by
            // *addr = 0;  /// important to write to this address so that a page would be mapped for the address by OS. this is _important_ for finding cha by
            // hashing method to not see the effects of lazy mapping. I DO NOT WANT TO WRITE IN THIS CASE!
        }
        else
        {
            //            SPDLOG_ERROR("addr is nullptr.");
            std::exit(EXIT_FAILURE);
        }

        const int cha = findCHAByHashing(reinterpret_cast<uintptr_t>(addr), base_sequence_28_skx);
        // const auto socket_cha = findCHAPerfCounter(reinterpret_cast<long long*>(addr));
        // std::cout << "i: " << i << std::endl;
        res[i] = cha;
        ++i;

        for (int j = 0; j < CACHE_LINE_SIZE / sizeof(T) - 1; ++j)
        {
            if (i < NUM)
            {
                // std::cout << "j: " << j << ", i: " << i << std::endl;
                // auto addr = static_cast<const int *>(&arr[i]);
                // const int cha_ = findCHAByHashing(reinterpret_cast<uintptr_t>(addr), base_sequence_18_skx);
                // assert(cha_ == cha);
                res[i] = cha;
                ++i;
            }
            else
            {
                break;
            }
        }

        //        SPDLOG_INFO("cha: {}", cha);
        // ++addr;
    }

    return res;
}

int main() {
    std::cerr << "here\n";
    if (geteuid() != 0) {
        std::cerr << "This program must be run as root!" << std::endl;
        return 1; // exit with an error code
    }
//#if 0
    std::cerr << "here1\n";
    std::cout << "Running with root privileges." << std::endl;
    std::cout << "NUM: " << NUM << std::endl;
    std::cout << "size in bytes: " << SIZE << std::endl;

    int destroy = 0;
    //printf("argv[1]: %s\n", argv[1]);
    //if(argc == 2 && strcmp("1", argv[1]) == 0) {
        //destroy = 1;
    //}
//#if 0
  // 4-byte shared memory 
  const int fd = shm_open(NAME/*"/dev/shm/shared_mem"*/, O_RDONLY, 0666);
  if (fd < 0) {
    perror("shm_open()");
    return EXIT_FAILURE;
  }
//#if 0
    char *const data = static_cast<char *const>(mmap(nullptr, SIZE, PROT_READ, MAP_SHARED, fd, 0));
    if (data == MAP_FAILED) {
        perror("mmap failed");
        std::exit(EXIT_FAILURE);
    }

  printf("consumer mapped address: %p\n", data);
//#if 0
#if 0
  for (int i = 0; i < NUM; ++i) {
    std::cout << "data[" << i << "]: " << data[i] << ", CHA of &data[" << i << "]: " << data[i + NUM] << std::endl; // read it here so that the address is mapped in so that physical address conversion does not mess up.
  }
#endif
  uint64_t holder = 0;
  for (int i = 0; i < NUM; ++i) {
    holder += (char) data[i];
  }
  std::cerr << holder << "\n";
//#if 0
  printf("physical address: %" PRIxPTR ", cha: %d\n", getPhysicalAddress(((uintptr_t)data)), findCHAByHashing(uintptr_t(data), base_sequence_28_skx));

    std::cout << "now benchmarking..." << std::endl;

    std::vector<char> shm_arr;
    long long elapsed_shm = -1;
    {
        auto t1 = high_resolution_clock::now();
        shm_arr = findWithShm<char>(reinterpret_cast<char*>(data));
        auto t2 = high_resolution_clock::now();
        auto ms_int = duration_cast<nanoseconds>(t2 - t1);
        elapsed_shm = ms_int.count();
        std::cout << "shm arr: " << elapsed_shm << " us\n";
    }

//#if 0
    std::vector<char> func_arr;
    long long elapsed_func = -1;
    {
        auto t1 = high_resolution_clock::now();
        func_arr = findWithHashFunc<char>(reinterpret_cast<char*>(data));
        auto t2 = high_resolution_clock::now();
        auto ms_int = duration_cast<nanoseconds>(t2 - t1);
        elapsed_func = ms_int.count();
        std::cout << "func arr: " << elapsed_func << " us\n";  
    }

#if 0
    std::array<int, NUM> func_all_arr;
    long long elapsed_func_all = -1;
    {
        auto t1 = high_resolution_clock::now();
        func_all_arr = findChasOfArray<int>(data);
        auto t2 = high_resolution_clock::now();
        auto ms_int = duration_cast<microseconds>(t2 - t1);
        elapsed_func_all = ms_int.count();
        std::cout << "func all arr: " << elapsed_func_all << " us\n";  
    }    
#endif
    std::cout << "diff: " << elapsed_func / static_cast<double>(elapsed_shm) << std::endl;
   // std::cout << "diff all: " << elapsed_func_all / static_cast<double>(elapsed_shm) << std::endl;

//#if 0
    std::cerr << "shm_arr.size(): " << shm_arr.size() << ", func_arr.size(): " << func_arr.size() << std::endl;
    int min_size = std::min(shm_arr.size(), func_arr.size());
    int i;
    for(i = 0; i < min_size; ++i) {
        if(shm_arr[i] != func_arr[i]) {
		std::cerr << "i: " << i << ", shm_arr[i]: " << (int) shm_arr[i] << ", func_arr[i]: " << (int) func_arr[i] << std::endl;
		break;
	}
        //assert(shm_arr[i] == func_all_arr[i]);
    }
    //std::cerr << "i: " << i << ", shm_arr[i]: " << (int) shm_arr[i] << ", func_arr[i]: " << (int) func_arr[i] << std::endl;
//#endif    

    if (i == min_size)
    	std::cout << "cha assert success." << std::endl;
    else
	std::cout << "fails.\n";
//#endif
    // 8-byte shared memory
#if 0
    const int fd1 = shm_open(NAME_DOUBLE, O_RDONLY, 0666);
  if (fd1 < 0) {
    perror("shm_open()");
    return EXIT_FAILURE;
  }

  double *const double_data = static_cast<double *const>(mmap(nullptr, DOUBLE_SIZE, PROT_READ, MAP_SHARED, fd1, 0));
    if (double_data == MAP_FAILED) {
        perror("mmap failed");
        std::exit(EXIT_FAILURE);
    }

  printf("consumer mapped address: %p\n", double_data);

  //uint64_t holder;
  for (int i = 0; i < NUM; ++i) {
    holder += double_data[i];
  }

  printf("physical address: %" PRIxPTR ", cha: %d, holder %0.2lf\n", getPhysicalAddress(((uintptr_t)double_data)), findCHAByHashing(uintptr_t(double_data), base_sequence_28_skx), holder);

    std::cout << "now benchmarking..." << std::endl;

    //std::array<double, NUM> shm_arr_double;
    long long elapsed_shm1 = -1;
    {   
        auto t1 = high_resolution_clock::now();
//        shm_arr = findWithShm<double>(double_data);
        auto t2 = high_resolution_clock::now();
        auto ms_int = duration_cast<microseconds>(t2 - t1);
        elapsed_shm1 = ms_int.count();
        std::cout << "shm arr: " << elapsed_shm1 << " us\n";
    }

    //std::array<int, NUM> func_arr;
    //long long elapsed_func = -1;
    {
        auto t1 = high_resolution_clock::now();
        func_arr = findWithHashFunc<double>(double_data);
        auto t2 = high_resolution_clock::now();
        auto ms_int = duration_cast<microseconds>(t2 - t1);
        elapsed_func = ms_int.count();
        std::cout << "func arr: " << elapsed_func << " us\n";
    }
#endif
    //std::cout << "diff: " << elapsed_func / static_cast<double>(elapsed_shm1) << std::endl;
   // std::cout << "diff all: " << elapsed_func_all / static_cast<double>(elapsed_shm) << std::endl;
#if 0
   for (int i = 0; i < 64; ++i) {
    std::cout << "shm_arr_double[" << i << "]: " << shm_arr_double[i] << ", func_arr[" << i << "]: " << func_arr[i] << std::endl;
  } 
#endif
#if 0
    for(int i = 0; i < NUM; ++i) {
        assert(shm_arr[i] == func_arr[i]);
        //assert(shm_arr[i] == func_all_arr[i]);
    }
#endif

    std::cout << "cha assert for double array success." << std::endl;
//#endif
  munmap(data, SIZE);
//#endif
  //munmap(double_data, DOUBLE_SIZE);
  close(fd);
  //close(fd1);
if(destroy) {
    printf("destroyed.\n");
    shm_unlink(NAME);
    //shm_unlink(NAME_DOUBLE);
}
//#endif  

  return EXIT_SUCCESS;
}

