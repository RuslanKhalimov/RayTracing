input_file_name = 'cmake-build-release/output.txt'
output_file_name = 'converted.radiance'

def process_table(text_file, out):
    line = text_file.readline()
    while line != '\n' and line != '':
        tokens = line.split()
        out.append(list(map(float, tokens)))
        line = text_file.readline()

data = [[], [], [], []]

with open(input_file_name) as text_file:
    line = text_file.readline()
    i = -1
    while line != '':
        tokens = line.split()
        if len(tokens) == 0:
            line = text_file.readline()
            continue
        if tokens[0] == 'wave_length':
            i += 1
        process_table(text_file, data[i])
        line = text_file.readline()

pp = PostProcessor(PPDataUnits.RADIANCE, [400, 500, 600, 700], *data)

pp.SaveToHDR(output_file_name, overwrite = OverwriteMode.OVERWRITE)

