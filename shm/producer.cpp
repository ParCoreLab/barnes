#include "protocol.h"

#include <cstdio>
#include <sys/mman.h>
#include <cinttypes>
#include <cstring>

#include <iostream>

// http://logan.tw/posts/2018/01/07/posix-shared-memory/
int main() {
    if (geteuid() != 0) {
        std::cerr << "This program must be run as root!" << std::endl;
        return 1; // exit with an error code
    }

    std::cout << "Running with root privileges." << std::endl;
    std::cout << "NUM: " << NUM << std::endl;
    std::cout << "SIZE in bytes: " << SIZE << std::endl;

// shared memory with 4-byte elements
  const int fd = shm_open(NAME, O_CREAT | O_EXCL | O_RDWR, 0600);
  if (fd < 0) {
    perror("shm_open()");
    return EXIT_FAILURE;
  }

  ftruncate(fd, SIZE);

  char *const data = static_cast<char *const>(mmap(nullptr, SIZE, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0));
    if (data == MAP_FAILED) {
        perror("mmap failed");
        std::exit(EXIT_FAILURE);
    }      

  for (int i = 0; i < NUM; ++i) {
    data[i] = 0;
    data[i + NUM] = findCHAByHashing((uintptr_t(&data[i])), base_sequence_28_skx);
    // std::cout << "data[" << i << "]: " << data[i] << ", CHA of &data[" << i << "]: " << data[i + NUM] << std::endl;
  }

  // print here so that physical address does not get printed as 0 as a result of the address not being mapped to the actual DRAM yet. by writing into the addresses
  // above, we make sure that mapping is done.
  printf("producer mapped virtual address: %p, physical address: %" PRIxPTR ", cha: %d\n", data,
    getPhysicalAddress(((uintptr_t)data)), findCHAByHashing((uintptr_t(data)), base_sequence_28_skx));

    // printf("waiting for char input...\n");
    // getchar();

  munmap(data, SIZE);
  close(fd);

#if 0
// shared memory with 8-byte elements
  //printf("before another shm_open\n"); 
  const int fd1 = shm_open(NAME_DOUBLE, O_CREAT | O_EXCL | O_RDWR, 0600);
  if (fd1 < 0) {
    perror("shm_open()");
    return EXIT_FAILURE;
  }

  //printf("before another ftruncate\n");
  ftruncate(fd1, DOUBLE_SIZE);

  //printf("before another mmap\n");  
  double *const double_data = static_cast<double *const>(mmap(nullptr, DOUBLE_SIZE, PROT_READ | PROT_WRITE, MAP_SHARED, fd1, 0));
   if (double_data == MAP_FAILED) {
        perror("mmap failed");
        std::exit(EXIT_FAILURE);
    }

  //printf("before another findCHAByHashing\n");
  for (int i = 0; i < NUM; ++i) {
    double_data[i] = 0;
    double_data[i + NUM] = (double) findCHAByHashing((uintptr_t(&double_data[i])), base_sequence_28_skx);
    //std::cout << "double_data[" << i << "]: " << double_data[i] << ", CHA of &double_data[" << i << "]: " << double_data[i + NUM] << std::endl;
  }

  //printf("after another findCHAByHashing\n");
  // print here so that physical address does not get printed as 0 as a result of the address not being mapped to the actual DRAM yet. by writing into the addresses
  // above, we make sure that mapping is done.
  printf("preducer mapped virtual address: %p, physical address: %" PRIxPTR ", cha: %d\n", double_data, 
    getPhysicalAddress(((uintptr_t)double_data)), findCHAByHashing((uintptr_t(double_data)), base_sequence_28_skx));

    // printf("waiting for char input...\n");
    // getchar();

  munmap(double_data, DOUBLE_SIZE);
  //printf("unmapped\n");
  close(fd1);
    //printf("closed fd1\n");
#endif

  return EXIT_SUCCESS;
}
