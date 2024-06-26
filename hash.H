//using namespace std;

class HashTable
{
	int INDEX_NUM;    // No. of buckets

	// Pointer to an array containing buckets
	std::pair<std::atomic<uint32_t>,std::list<std::pair<long unsigned int, std::pair<int, int>>>> *table;
	public:
	HashTable(int index_num);  // Constructor
	HashTable();
	~HashTable();  // Destructor

	// inserts a key into hash table
	void insertItem(long unsigned int index, long unsigned int cacheLine);

	bool findItem(long unsigned int index, long unsigned int cacheLine);

	void wait(long unsigned int index);

	void signal(long unsigned int index);

	std::pair<long unsigned int, std::pair<int, int>>& getItem(long unsigned int index, long unsigned int cacheLine);

	// hash function to map values to key
	long unsigned int hashFunction(long unsigned int cacheLine) {
		return (cacheLine % INDEX_NUM);
	}

	void displayHash();

	void print();
	void print_long(long unsigned int index);

	HashTable& operator=(HashTable other)
    	{
        	std::swap(INDEX_NUM, other.INDEX_NUM);
        	std::swap(table, other.table);
        	return *this;
    	}
};

HashTable::HashTable()
{
        INDEX_NUM = 0;
}

HashTable::HashTable(int index_num)
{
	std::cerr << "in constructor begins\n";
	INDEX_NUM = index_num;
	table = new std::pair<std::atomic<uint32_t>,std::list<std::pair<long unsigned int, std::pair<int, int>>>>[INDEX_NUM];
	for(int i = 0; i < INDEX_NUM; i++)
		table[i].first = 1;
	//for(int i = 0; i < INDEX_NUM; i++)
          //      std::cerr << table[i].first << "\n";
	std::cerr << "in constructor ends\n";
}

HashTable::~HashTable()
{
	delete[] table;
}

void HashTable::insertItem(long unsigned int index, long unsigned int cacheLine)
{
	//long unsigned int index = hashFunction(cacheLine);
	table[index].second.push_back({cacheLine, {0, 0}}); 
}

#if 0
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
#endif
bool HashTable::findItem(long unsigned int index, long unsigned int cacheLine)
{
	// get the hash index of key
	//long unsigned int index = hashFunction(cacheLine);

	// find the key in (index)th list
	std::list <std::pair<long unsigned int, std::pair<int, int>>> :: iterator i;
	for (i = table[index].second.begin();
			i != table[index].second.end(); i++) {
		if (i->first == cacheLine)
			break;
	}

	// if key is found in hash table, remove it
	if (i != table[index].second.end())
		return true;
	return false;
}

void HashTable::wait(long unsigned int index)
{
        //while(table[index].first <= 0);
	//table[index].first--;
	auto oldval = table[index].first.load();
        while (oldval == 0 || !table[index].first.compare_exchange_strong(oldval, 0)) {
           oldval = table[index].first.load();
	   //std::cerr << "in wait, val: " << oldval << "\n";
	}
}

void HashTable::signal(long unsigned int index)
{
	//std::cerr << "in signal\n";
        auto val = table[index].first.fetch_add(1, std::memory_order_seq_cst);
	//std::cerr << "in signal, val: " << val << "\n";
}

void HashTable::print_long(long unsigned int index)
{
	std::cerr << "hello world long\n";
}

void HashTable::print()
{
        std::cerr << "hello world\n";
        //auto val = table[index].first.fetch_add(1, std::memory_order_seq_cst);
        //std::cerr << "in signal, val: " << val << "\n";
}

std::pair<long unsigned int, std::pair<int, int>>& HashTable::getItem(long unsigned int index, long unsigned int cacheLine)
{
	// get the hash index of key
	//long unsigned int index = hashFunction(cacheLine);

	// find the key in (index)th list
	std::list <std::pair<long unsigned int, std::pair<int, int>>> :: iterator i;
	for (i = table[index].second.begin();
			i != table[index].second.end(); i++) {
		if (i->first == cacheLine)
			break;
	}

	// if key is found in hash table, remove it
	if (i != table[index].second.end()) {
		return *i;//table[index];
	}
	std::pair<long unsigned int, std::pair<int, int>> empty_pair({0, {0, 0}});
	return empty_pair;
}

// function to display hash table
void HashTable::displayHash() {
	for (int i = 0; i < INDEX_NUM; i++) {
		std::cout << i;
		for (auto x : table[i].second)
			std::cout << " --> cache line: " << x.first << ", a: " << x.second.first << ", b: " << x.second.second;
		std::cout << std::endl;
	}
}
