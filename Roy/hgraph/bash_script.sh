ls *|sed -r 's/^tensors-(.*).pkl$/\1/'|awk '$1>40 {print "tensors-"$1".pkl"}'| xargs rm