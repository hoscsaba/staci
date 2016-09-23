import matplotlib.pyplot as plt
import sys

infile = open(str(sys.argv[1]))
print("Datafile plotted: %s" % str(sys.argv[1]))

vec_fun_eval = []
vec_obj_val = []
fun_eval_prev = 0
for line in infile:
    parts = line.split(';')
    fun_eval = int(parts[0])
    obj_val = float(parts[1])
    vec_fun_eval.append(fun_eval_prev)
    vec_obj_val.append(obj_val)
    vec_fun_eval.append(fun_eval)
    vec_obj_val.append(obj_val)
    fun_eval_prev = fun_eval

plt.plot(vec_fun_eval, vec_obj_val)
plt.xlabel("function evaluation")
plt.ylabel("objective")
plt.grid(True)

plt.show()
#plt.savefig("tmp.pdf")
