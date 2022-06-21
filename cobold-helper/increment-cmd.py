#!/home/sdeshmukh/anaconda3/bin/python
# Increment the 'cobold.cmd' file by a value passed in, writing to a new
# 'cobold-new.cmd'
with open('./cobold.cmd', 'r', encoding='utf-8') as infile:
  lines = infile.readlines()

split_line = lines[-1].split('.')
model_string, identifier = split_line[:-1], split_line[-1].split(' ')
prev_final_val, filler = int(identifier[0]), ' '.join(
    identifier[1:]).replace('*', ' ')
start_val = prev_final_val + 1
increment_value = 40  # Increment each entry by this

new_lines = []
for i in range(start_val, start_val + increment_value):
  new_value = str(i).zfill(4)
  new_lines.append(''.join(model_string) + f".{new_value} {filler}")

with open('cobold-new.cmd', 'w', encoding='utf-8') as outfile:
  outfile.writelines(new_lines)
