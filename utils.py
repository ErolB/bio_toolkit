def read_fasta(path):
    # labels are the header text after '>', sequences are joined across wrapped lines
    sequences = {}
    label = None
    chunks = []
    with open(path, encoding='utf-8') as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.strip()
            if not line or line.startswith(';'):  # blank line or legacy comment
                continue
            if line.startswith('>'):
                if label is not None:
                    sequences[label] = ''.join(chunks)
                label = line[1:].strip()
                if not label:
                    raise ValueError(f'{path} line {line_number}: header has no label')
                if label in sequences:
                    raise ValueError(f'{path} line {line_number}: duplicate label {label!r}')
                chunks = []
            else:
                if label is None:
                    raise ValueError(f'{path} line {line_number}: sequence data before any header')
                chunks.append(''.join(line.split()))  # drop any internal whitespace
    if label is not None:
        sequences[label] = ''.join(chunks)
    return sequences
