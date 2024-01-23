Traceback (most recent call last):
  File "./train_translator.py", line 52, in <module>
    vocab = [x.strip("\r\n ").split() for x in open(args.vocab)]
  File "./train_translator.py", line 52, in <listcomp>
    vocab = [x.strip("\r\n ").split() for x in open(args.vocab)]
OSError: [Errno 5] Input/output error
