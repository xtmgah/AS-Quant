class C:
	def __init__(self, value, index):
		self.value = value
		self.index = index

	def __cmp__(self, other):
		return cmp(self.value, other[self.index])

class EXON:
	def __init__(self):
		self.st = 0
		self.en = 0

class Stack:
	def __init__(self):
		self.items = []

	def size(self):
		return len(self.items)

	def isEmpty(self):
		return self.items == []

	def push(self, val):
		self.items.append(val)

	def top(self):
		if self.isEmpty():
			return None
		else:
			return self.items[self.size()-1]

	def pop(self):
		if self.isEmpty():
			return None
		else:
			return self.items.pop()