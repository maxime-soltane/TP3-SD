import hashlib

class SimpleBloomFilter:
	def __init__(self, size=100, num_hashes=1):
		self.size = size
		self.num_hashes = num_hashes
		self.bit_array = [0] * size

	def _hashes(self, item):
		hash_values = []
		for i in range(self.num_hashes):
			hash_func = hashlib.sha256((str(i) + item).encode()).hexdigest()
			hash_values.append(int(hash_func, 16) % self.size)
		return hash_values

	def add(self, item):
		for pos in self._hashes(item):
			self.bit_array[pos] = 1

	def contains(self, item):
		return all(self.bit_array[pos] for pos in self._hashes(item))

	def merge(self, other):
		assert self.size == other.size, "Bloom filters must be of the same size!"
		merged_filter = SimpleBloomFilter(self.size, self.num_hashes)
		merged_filter.bit_array = [a | b for a, b in zip(self.bit_array, other.bit_array)]
		return merged_filter


class StructureNode:
	def __init__(self, bloom_filter=None):
		self.bloom = bloom_filter if bloom_filter else SimpleBloomFilter()
		self.left = None
		self.right = None
		self.datasets = []  # list of dataset names at leaf nodes

class Structure:
	def __init__(self, datasets, kmers_dict, bloom_size=10000, num_hashes=3):
		self.leaves = {}  # maps dataset names to their Bloom filter nodes
		self.root = self._build_tree(datasets, kmers_dict, bloom_size, num_hashes)

	def _build_tree(self, datasets, kmers_dict, bloom_size, num_hashes):
		nodes = []

		# Step 1
		for dataset in datasets:
			bf = SimpleBloomFilter(bloom_size, num_hashes)
			for kmer in kmers_dict[dataset]:
				bf.add(kmer)
			node = StructureNode(bf)
			node.datasets = [dataset]
			self.leaves[dataset] = node
			nodes.append(node)

		# Step 2
		while len(nodes) > 1:
			new_nodes = []
			for i in range(0, len(nodes), 2):
				if i + 1 < len(nodes): 
					merged_bf = nodes[i].bloom.merge(nodes[i + 1].bloom)
					parent = StructureNode(merged_bf)
					parent.left = nodes[i]
					parent.right = nodes[i + 1]
					parent.datasets = nodes[i].datasets + nodes[i + 1].datasets
				else:
					parent = nodes[i] 
				new_nodes.append(parent)
			nodes = new_nodes

		return nodes[0] if nodes else None  

	def query(self, kmer):
		results = []
		self._query_recursive(self.root, kmer, results)
		return results

	def _query_recursive(self, node, kmer, results):
		if node is None:
			return
		if node.bloom.contains(kmer): 
			if node.left is None and node.right is None: 
				results.extend(node.datasets)
			else:
				self._query_recursive(node.left, kmer, results)
				self._query_recursive(node.right, kmer, results)


datasets = ["Dataset1", "Dataset2", "Dataset3", "Dataset4"]
kmers_dict = {
	"Dataset1": ["ACGT", "TGCA", "GCTA"],
	"Dataset2": ["CGTA", "GCTA", "TACC"],
	"Dataset3": ["AAGT", "TCCA", "CGGT"],
	"Dataset4": ["TGGC", "GGCA", "CCAA"]
}
#test

structure = Structure(datasets, kmers_dict, bloom_size=100, num_hashes=1)
query_kmers = ["GCTA", "TCCA", "ACGT", "GGGG"]
for kmer in query_kmers:
	result = structure.query(kmer)
	print(f"K-mer '{kmer}' found in datasets: {result}")
