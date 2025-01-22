from typing import List, Dict, Optional, Tuple
import re

SUMMARY: str = ''

def determine_sequence_type(sequence: str) -> str:
    # using negative look ahead (?!) to ensure the character is not present at any position of the sequence.
    # NOTE: referred to the wikipedia page on FASTA (https://en.wikipedia.org/wiki/FASTA_format) for determining the ambiguous characters
    valid_dna: str = r'(?!.*U)[ACGTNRYKMSWBDHVXZ]+'
    valid_rna: str = r'(?!.*T)[ACGUNRYKMSWBDHVXZ]+'
    
    # Determine sequence type
    # Ambiguous nucleotides for DNA/RNA = [NRYKMSWBDHVX]
    if re.match(valid_dna, sequence):
        return 'DNA'
    elif re.match(valid_rna, sequence):
        return 'RNA'
    else:
        return 'Invalid'

def identify_orf(sequence: str) -> Tuple[str, int]:
    count = 0
    sequence_type = determine_sequence_type(sequence)

    # Unsure BUG: is it necessary for the start codons to start in the multiples of 3? The example provided says otherwise.
    # In this implementation, I am assuming that the start codon can begin at any index, but thereafter the checking happens STRICTLY in steps of 3.
    # So the stop codon will have to begin at one of the indices [startIndex + 3*r], where r = [1, 2, 3, ...]
    # 
    # 
    # If we go by this article on Wikipedia (https://en.wikipedia.org/wiki/Open_reading_frame), the sequence is divided
    # into characters of 3 and the possible indices for start codon is 0, 3, 6, etc. If we want the implementation to go
    # by this it just requires a if statememt to check if startIndex % 3 == 0? If yes, we accept, else reject.

    orfs = ''
    start = 0

    if sequence_type == 'DNA':
        stop_codons = ['TAA', 'TAG', 'TGA']
        start_codon = 'ATG'
    elif sequence_type == 'RNA':
        stop_codons = ['UAA', 'UAG', 'UGA']
        start_codon = 'AUG'
    else:
        ans = f'Found {count} ORF:\n' + orfs
        return (ans, count)

    while start < len(sequence):

            # Look for the first 'ATG' from start index 
            matchIndex = sequence.find(start_codon, start)
            
            # If ATG not found, break loop
            if matchIndex == -1:
                break

            # If match found, from that index keep checking every 3 steps for stop codons
            for r in range(1, ((len(sequence) - matchIndex)// 3) + 1):
                stop_codon_starting_index = matchIndex + r*3
                stop_codon_ending_index = stop_codon_starting_index + 3
                # print(f'\n\nStop index: {stop_codon_ending_index}')

                if stop_codon_ending_index > len(sequence):
                    break

                # print(f'Checking substring: {sequence[matchIndex: stop_codon_ending_index]}')

                if sequence[stop_codon_starting_index: stop_codon_ending_index] in stop_codons:
                    # print(f'Valid: {sequence[matchIndex: stop_codon_ending_index]}')
                    count += 1
                    orfs += f'- Start at {matchIndex}, Stop at {stop_codon_starting_index}, ORF: {sequence[matchIndex: stop_codon_ending_index]}\n'
                
            start = matchIndex + 3

        
    ans = f'Found {count} ORF:\n' + orfs

    return (ans, count)


def compute_sequence_stats(sequence: str, sequence_type: str) -> str:
    length: str = len(sequence)
    counts: Dict[str, int] = {'A': 0, 'C': 0, 'G': 0}
    (_, orf_count) = identify_orf(sequence)

    if sequence_type == 'DNA':
        counts['T'] = 0
    elif sequence_type == 'RNA':
        counts['U'] = 0
    else:
        return ''
    
    counts['Ambiguous'] = 0

    for char in sequence:
        if char in ['A', 'C', 'T', 'G', 'U']:
            counts[char] += 1
        else:
            counts['Ambiguous'] += 1
    
    return f'Length: {length}\nCounts: {"".join([f"{key}: {val}," for key, val in counts.items()])}\nNumber of ORFs: {orf_count}\n'



def extract_sequence(group: str) -> str:
    if not group: 
        raise ValueError('group is of `NoneType`')
    
    # Will use regex and extract the sequence using groups, ie. after the header line
    sequence: str = re.search(r'>(.+)\n([A-Za-z]+)\*?', group).group(2)
    sequence = sequence.upper()
    return sequence


def find_type(group: str) -> str:
    sequence = extract_sequence(group)
    return determine_sequence_type(sequence)
    

def fasta_parser(filename: Optional[str]) -> List[str] | None:
    """Performs the deliverable for Part 1"""

    if not filename:
        raise FileNotFoundError(f'Could not find file: {filename}')
    
    groups = []
    total_dna_length = 0
    dna_count = 0
    total_rna_length = 0
    rna_count = 0
    invalid_count = 0
    global SUMMARY

    try:
        with open(filename, 'r') as file:
            group: str = ''
            for line in file:
                line: str = line.strip(' \n\t')
                res = re.search('^>', line)

                if not res:
                    group += line
                else:
                    if group:
                        sequence = extract_sequence(group)
                        seq_type = find_type(group)
                        group += f'\nType: {seq_type}\n'

                        if filename == 'test_files/test_part2_in.fasta':
                            (orfs, _) = identify_orf(sequence)
                            group +=  orfs 
                        elif filename == 'test_files/test_part3_in.fasta':
                            stats = compute_sequence_stats(sequence, seq_type)
                            group += stats

                            if seq_type == 'DNA':
                                dna_count += 1
                                total_dna_length += len(sequence)
                            elif seq_type == 'RNA':
                                rna_count += 1
                                total_rna_length += len(sequence)
                            else:
                                invalid_count += 1
                        
                        groups.append(group)

                    group = '' + line + '\n'

            if filename == 'test_files/test_part3_in.fasta':
                SUMMARY = f"""\n---SUMMARY---\nValid DNA sequences: {dna_count}\nMean DNA length: {total_dna_length/dna_count}\nValid RNA sequences: {rna_count}\nMean RNA Length: {total_rna_length/rna_count}\nInvalid Sequences: {invalid_count}\n"""
            
            return groups
        
    except Exception as e:
        raise RuntimeError(e)

def main():
    """Main function to run the FASTA parser and display output"""
    try:
        print('--- OUTPUT for PART 1 ----\n')
        groups = fasta_parser('test_files/test_part1_in.fasta')
        for group in groups:
            print(group)

        print('\n\n--- OUTPUT for PART 2 ----\n')
        groups = fasta_parser('test_files/test_part2_in.fasta')
        for group in groups:
            print(group)

        print('\n\n--- OUTPUT for PART 3 ----\n')
        groups = fasta_parser('test_files/test_part3_in.fasta')
        for group in groups:
            print(group)
        print(SUMMARY)

    except Exception as e:
        print(f'An error occurred: {e}')

if __name__ == "__main__":
    main()