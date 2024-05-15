## Requirements

* A 28-core per socket Intel Cascade Lake machine running a Linux OS kernel

## How to prepare the CHA-aware memory manager

```bash
$ cd shm
$ make
$ sudo ./automate.sh
```

## How to compile

```bash
$ make
```

## How to run the CHA-aware BARNES on an input configuration

```bash
$ sudo ./BARNES < input
```

## How to run the traffic comparison experiment

```bash
sudo ./input_size_compare.sh
```
