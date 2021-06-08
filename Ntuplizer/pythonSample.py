import sys

print('Number of arguments:', len(sys.argv), 'arguments.')
print('Argument List:', str(sys.argv))
print('Second argument:', float(sys.argv[2]))

nLoops = int(sys.argv[2])

for i in range(nLoops):
    text = str(sys.argv[1])
    print(text)
