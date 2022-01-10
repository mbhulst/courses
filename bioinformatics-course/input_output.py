# TO LOAD INPUT:
with open('data.txt', 'r') as data:  # read data into string
    x = data.read()

with open('data.txt', 'r') as data:  # if the data is split by enters
    x = data.read().splitlines()

with open('data.txt', 'r') as data:  # if the data is split by spaces
    x = data.read().split()

with open('data.txt', 'r') as data:  # if the data is split by enters and should be converted to floats
    data = data.read().splitlines()
    data = [x.split() for x in data]
    data = [[float(x) for x in y] for y in data]

with open('data.txt', 'r') as file:  # if data is given as graph
    graph = dict(line.strip().split(' -> ') for line in file)
for key in graph:
    graph[key] = graph[key].split(',')
graph = {int(k): [int(i) for i in v] for k, v in graph.items()}


# TO FORMAT ANSWER
print(' '.join(map(str, answer)))   # if the result needs to be split by spaces

print('\n'.join(map(str, answer)))  # if the result needs to be split by enters

print('->'.join(map(str, answer)))  # if the result needs to be split by an arrow

for key in sorted(answer):
    print(key + " -> " + ",".join(answer[key]))
