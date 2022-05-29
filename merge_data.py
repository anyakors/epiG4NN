'''
Make input from intersection data
'''


import argparse
from tqdm import tqdm


def check_argv():

    parser = argparse.ArgumentParser(description=__doc__.strip().split('\n')[0], add_help=False)
    parser.add_argument('-m', '--method', help='method')
    parser.add_argument('--read', help='read filename')
    parser.add_argument('--write', help='write filename')

    return parser.parse_args()


def scores(read_fn, write_fn):

    def write_entry(fw, *args):

        avg_score = str(round(sum(args[-1]) / len(args[-1]), 6))
        args = list(args)
        args[-1] = avg_score

        fw.write('\t'.join(args) + '\n')

    with open(read_fn, 'r') as fr, open(write_fn, 'w') as fw:

        curr_id, curr_strand, curr_scores = None, None, []

        for record in tqdm(iter(lambda: fr.readline(), '')):

            if len(record.rstrip('\n').split('\t'))==10:
                pqs_chr, pqs_start, pqs_end, pqs_id, pqs_none, pqs_strand, exp_chr, exp_start, exp_end, exp_score = record.rstrip('\n').split('\t')
            else:
                pqs_chr, pqs_start, pqs_end, pqs_id, pqs_strand, exp_chr, exp_start, exp_end, exp_score = record.rstrip('\n').split('\t')
            pqs_start, pqs_end, exp_start, exp_end = map(int, (pqs_start, pqs_end, exp_start, exp_end))
            exp_score = float(exp_score)

            assert pqs_chr == exp_chr and exp_start < pqs_end and pqs_start < exp_end

            if pqs_id != curr_id or pqs_strand != curr_strand:
                if curr_id is not None:
                    write_entry(fw, curr_id, curr_strand, curr_scores)
                curr_id, curr_strand, curr_scores = pqs_id, pqs_strand, []

            curr_scores.append(exp_score)

        write_entry(fw, curr_id, curr_strand, curr_scores)


def meth_faire(read_fn, write_fn):

    with open(read_fn, 'r') as fr, open(write_fn, 'w') as fw:

        curr_pos = None

        for record in tqdm(iter(lambda: fr.readline(), '')):

            if len(record.rstrip('\n').split('\t'))==10:
                pqs_chr, pqs_start, pqs_end, pqs_id, pqs_none, pqs_strand, exp_chr, exp_start, exp_end, exp_score = record.rstrip('\n').split('\t')
            else:
                pqs_chr, pqs_start, pqs_end, pqs_id, pqs_strand, exp_chr, exp_start, exp_end, exp_score = record.rstrip('\n').split('\t')
            pqs_start, pqs_end, exp_start, exp_end = map(int, (pqs_start, pqs_end, exp_start, exp_end))

            assert pqs_chr == exp_chr and exp_start < pqs_end and pqs_start < exp_end

            fw.write(
                '\t'.join((
                    pqs_chr, 
                    str(curr_pos if curr_pos is not None else pqs_start), 
                    str(min(pqs_end, exp_end)), 
                    pqs_id, 
                    pqs_strand, 
                    exp_score
                )) + '\n'
            )
            
            curr_pos = exp_end if exp_end < pqs_end else None


if __name__ == '__main__':

    import os
    cwd = os.path.dirname(os.path.realpath(__file__))
    data_dir = os.path.join(cwd)#, os.pardir, 'data_new/293T'

    args = check_argv()
    read_fn = os.path.join(data_dir, args.read)
    write_fn = os.path.join(data_dir, args.write)

    if args.method == 'scores':
        scores(read_fn, write_fn)
    elif args.method == 'meth' or args.method == 'faire':
        meth_faire(read_fn, write_fn)
    else:
        print('Unidentified method argument. Exiting.')