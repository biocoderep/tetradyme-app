import json
from core import process_primers

def main():
    with open('sequences.fasta', 'r', encoding='utf-8') as f:
        fasta_content = f.read()
    with open('example.csv', 'r', encoding='utf-8') as f:
        csv_content = f.read()

    params = {
        'na': 50.0,
        'mg': 0.0,
        'dntps': 0.0,
        'tm_tolerance': 1.0,
        'of_primer': '',
        'or_primer': ''
    }

    result = process_primers(fasta_content, csv_content, params)
    print("Success Count:", result['success_count'])
    print("Failed Count:", len(result['failed_snps']))
    print("First result:", json.dumps(result['results'][0] if result['results'] else {}, indent=2))

if __name__ == '__main__':
    main()
