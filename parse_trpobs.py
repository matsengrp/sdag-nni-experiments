import sys

if (len(sys.argv) != 2):
    print("Usage: <trprobs_file>")
    exit(0)
else:
    filename = sys.argv[1]

fp = open(filename, "r")

for line in fp:
    # print("line: ", line)
    words = line.split()
    if (len(words) > 0 and words[0] == "tree"):
        print(words[len(words) - 1])

fp.close()
