import sys
import util2
from subprocess import call

state, state_int = util2.read_state(sys.argv[1], to_elements=False)
particle_ids = state_int[state_int[:, 1] == 0, 0]

batch_size = 256
call(["mkdir", "tmp", "-p"])

def do_work(i, start_index, length):
	ids = particle_ids[start_index:start_index+length]
	print(" ".join([str(x) for x in ids]))
	call(["bin/prune-track", sys.argv[2], "tmp/pruned", "-w{0}".format(",".join([str(x) for x in ids]))])
	call(["bin/export-track", "--precision", "5", "tmp/pruned", "tmp/exported"])
	call(["script/resamp.pl", "--if", "tmp/exported", "--of", "tmp/out{0}".format(i), "--tp", "{0}".format(particle_ids.max() + 1), "--p", "5", "--q", "2", "--pl", "4"])

for i in range(len(particle_ids) // batch_size):
	do_work(i, i * batch_size, batch_size)
	finalindex = (i + 1) * batch_size
	n = i

if finalindex != len(particle_ids):
	n = i + 1
	do_work(n, finalindex, len(particle_ids) - finalindex)

with open(sys.argv[3], "w") as out:
	for i in range(n+1):
		with open("tmp/out{0}".format(i), "r") as f:
			for line in f:
				if len(line) > 0 and line[0] == "#":
					continue
				out.write("{0}".format(line))

