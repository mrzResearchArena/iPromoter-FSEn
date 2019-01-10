def FASTA(filename):
    try:
        f = open(filename)
    except IOError:
        print ("The file, %s, does not exist" % filename)
        return

    order = []
    sequences = {}

    for line in f:
        if line.startswith('>'):
            name = line[1:].rstrip('\n')
            name = name.replace('_',' ')
            order.append(name)
            sequences[name] = ''

        else:
            sequences[name] += line.rstrip('\n').rstrip('*')

    print ("%d sequences found" % len(order))
    return order, sequences