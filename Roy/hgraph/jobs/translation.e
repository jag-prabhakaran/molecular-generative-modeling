Traceback (most recent call last):
  File "./translate.py", line 46, in <module>
    vocab = [x.strip("\r\n ").split() for x in open(args.vocab)]
  File "./translate.py", line 46, in <listcomp>
    vocab = [x.strip("\r\n ").split() for x in open(args.vocab)]
OSError: [Errno 5] Input/output error
