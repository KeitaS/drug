#coding:utf-8
"""
riboモデルに入力したDoseを調節するための係数を調べたもの
"""
import re, os

fname = "result/ribo.csv"
max = [0, 0]
min = [0, 0]
data = []

with open(fname, "r") as rf:
    for index, line in enumerate(rf.readlines()):
        line = re.split(",", line.strip())
        line = map(float, line)
        if index == 0:
            max = line
            min = line
        else:
            if max[1] < line[1]:
                max = line
            if min[1] > line[1]:
                min = line

        data.append(line)

ave = (max[1] + min[1]) / 2
# print max, min, ave

sa = 0
result = []
for index, v in enumerate(data):
    if index == 0 or abs(ave - v[1]) < sa:
        sa = abs(ave - v[1])
        result = v

print "result: [%e, %e]" % (result[0], result[1])
print "max: [%e, %e]" % (max[0], max[1])
print "min: [%e, %e]" % (min[0], min[1])
print "averate: %e" % ave
print "sa: %e" % sa

