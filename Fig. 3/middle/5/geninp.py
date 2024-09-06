# Read the original file
with open('ori.txt', 'r') as file:
    lines = file.readlines()
    num_lines = len(lines)
    print(num_lines)

# Extract the first five columns and create 'real.txt'
with open('real.txt', 'w') as real_file:
    for line in lines:
        columns = line.strip().split(',')
        first_five_columns = ','.join(columns[:9])
        real_file.write(first_five_columns + '\n')

# Get the number of lines
num_lines = len(lines)
num_lines = 5145-31+1
# Add the number of lines to the first column and create 'fake.txt'
with open('fake.txt', 'w') as fake_file:
    for i, line in enumerate(lines):
        columns = line.strip().split(',')
        columns[0] = str(int(columns[0]) + num_lines)
        fake_file.write(','.join(columns) + '\n')
