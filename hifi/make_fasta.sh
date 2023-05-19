module load seqtk/1.3

awk '/^S/{print ">"$2"\n"$3}' $1 | seqtk seq -l 80
