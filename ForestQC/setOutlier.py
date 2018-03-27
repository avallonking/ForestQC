import gzip
import sys
import os
import linecache

# def set_outlier(file_list):
#   DP_GQ_START_IDX = 9
#   DP_IDX = 2
#   GQ_IDX = 3
#   DP_THRESHOLD = 0
#   GQ_THRESHOLD = 0
#   gq = []
#   dp = []
#   # sample_idx = None
#   for file in file_list:
#     f = gzip.open(file, 'rt') if 'gz' in file else open(file, 'r')
#     for line in f:
#       if not line.startswith('#'):
#         variant_info = line
#         individual_info = variant_info.split('\t')[DP_GQ_START_IDX:]
#         # if not sample_idx:
#           # sample_idx = random.sample(range(len(individual_info)), int(len(individual_info) / 10))
#         # individual_info = itemgetter(*sample_idx)(individual_info)
#         for idv in individual_info:
#           idv_info = idv.split(':')
#           if len(idv_info) < 4:
#             continue
#           if idv_info[DP_IDX] != '.':
#             dp.append(float(idv_info[DP_IDX]))
#           if idv_info[GQ_IDX] != '.':
#             gq.append(float(idv_info[GQ_IDX]))
#
#     f.close()
#
#   DP_THRESHOLD = np.percentile(dp, 25, interpolation = 'lower')
#   GQ_THRESHOLD = np.percentile(gq, 25, interpolation = 'lower')
#   print('Outlier_DP: ' + str(DP_THRESHOLD))
#   print('Outlier_GQ: ' + str(GQ_THRESHOLD))


class VCFProcessor(object):
    BLOCK_FILENAME_FORMAT = '{0}_block_{1}.dat'

    def __init__(self, filenames):
        # filenames should be a list of vcf files
        self.filenames = filenames
        self.gq_block_filenames = []
        self.dp_block_filenames = []

    def write_block(self, data, _type, block_number):
        data.sort()
        filename = self.BLOCK_FILENAME_FORMAT.format(_type, block_number)
        file = open(filename, 'w')
        file.write('\n'.join(list(map(str, data))))
        file.close()
        if _type == 'gq':
            self.gq_block_filenames.append(filename)
        elif _type == 'dp':
            self.dp_block_filenames.append(filename)

    def get_block_filenames(self):
        return self.gq_block_filenames, self.dp_block_filenames

    def get_gq_dp(self, block_size):
        gq = []
        dp = []
        gq_block_num = 0
        dp_block_num = 0

        for filename in self.filenames:
            if '.gz' in filename:
                f = gzip.open(filename, 'rt')
            else:
                f = open(filename, 'r')

            for line in f:
                if not line.startswith('#'):
                    info = line.split('\t')[9:]
                    for sample_info in info:
                        sample_info = sample_info.split(':')
                        if len(sample_info) < 4:
                            continue
                        try:
                            gq.append(int(sample_info[3]))
                        except ValueError:
                            pass
                        try:
                            dp.append(int(sample_info[2]))
                        except ValueError:
                            pass

                if sys.getsizeof(gq) >= block_size:
                    self.write_block(gq, 'gq', gq_block_num)
                    gq = []
                    gq_block_num += 1
                if sys.getsizeof(dp) >= block_size:
                    self.write_block(dp, 'dp', dp_block_num)
                    dp = []
                    dp_block_num += 1
            f.close()

        self.write_block(gq, 'gq', gq_block_num)
        self.write_block(dp, 'dp', dp_block_num)

        return gq_block_num, dp_block_num

    def cleanup(self):
        any(map(lambda f: os.remove(f), self.gq_block_filenames))
        any(map(lambda f: os.remove(f), self.dp_block_filenames))


class NWayMerge(object):
    def select(self, choices):
        min_index = -1
        # while min_index not in choices:
        # find the smallest index that still in choices.
        # But remember, choice is automatically sorted by string in a dict,
        # so the order is: 0,1,10,11,12,13,14,2,3,4,5,6,7,8,9
        #    min_index += 1

        # ensure min_index is in choices
        for an_index in choices:
            min_index = an_index
            break

        for i in choices:
            if int(choices[i].strip()) < int(choices[min_index].strip()):
                min_index = i

        return min_index


class FilesArray(object):
    def __init__(self, files):
        self.files = files
        self.empty = set()
        self.num_buffers = len(files)
        self.buffers = {i: None for i in range(self.num_buffers)}

    def get_dict(self):
        return {i: self.buffers[i] for i in range(self.num_buffers) if i not in self.empty}

    def refresh(self):
        for i in range(self.num_buffers):
            if self.buffers[i] is None and i not in self.empty:
                self.buffers[i] = self.files[i].readline()

                if self.buffers[i] == '':
                    self.empty.add(i)

        if len(self.empty) == self.num_buffers:
            return False

        return True

    def unshift(self, index):
        value = self.buffers[index]
        self.buffers[index] = None
        return value.strip()


class FileMerger(object):
    def __init__(self, merge_strategy):
        self.merge_strategy = merge_strategy

    def get_file_handles(self, filenames, buffer_size):
        files = {}
        for i in range(len(filenames)):
            files[i] = open(filenames[i], 'r', buffer_size)
        return files

    def close_files(self, file_handles):
        map(lambda f: f.close(), file_handles)

    def merge(self, filenames, outfilename, buffer_size):
        outfile = open(outfilename, 'w', buffer_size)
        file_handles = self.get_file_handles(filenames, buffer_size)
        buffers = FilesArray(file_handles)
        while buffers.refresh():
            min_index = self.merge_strategy.select(buffers.get_dict())
            outfile.write(buffers.unshift(min_index) + '\n')
        self.close_files(file_handles)
        outfile.close()


class ExternalSort(object):
    def __init__(self, block_size):
        self.block_size = block_size

    # def get_number_blocks(self, filenames, block_size):
    #     total_size = 0
    #     for filename in filenames:
    #         total_size += os.stat(filename).st_size
    #     return int(total_size / block_size) + 1

    def sort(self, filenames, outfilename):
        # filenames is a list of vcf file(s)
        vcf_processor = VCFProcessor(filenames)
        gq_block_num, dp_block_num = vcf_processor.get_gq_dp(self.block_size)
        gq_block_filenames, dp_block_filenames = vcf_processor.get_block_filenames()

        merger = FileMerger(NWayMerge())
        # assuming we only have memory size of block_size
        # share the memory with all files needed to be open (block files and outfile)
        gq_buffer_size = int(self.block_size / (gq_block_num + 1))
        dp_buffer_size = int(self.block_size / (dp_block_num + 1))

        # merge gq and dp separately
        merger.merge(gq_block_filenames, 'gq.' + outfilename, gq_buffer_size)
        merger.merge(dp_block_filenames, 'dp.' + outfilename, dp_buffer_size)

        vcf_processor.cleanup()


# use external merge sort to sort the quality vector first, then take quartile out
def get_quartile(filename):
    total_line_count = sum(1 for line in open(filename, 'r'))
    # the line of quartile number
    quartile_line = int(int(total_line_count) * 0.25)
    return linecache.getline(filename, quartile_line).strip()


def set_outlier(filenames, temp_dir, outfilename, block_size):
    # filenames should be a list of vcf files
    if block_size[-1].lower() == 'g':
        block_size = int(block_size[:-1]) * 1024 * 1024 * 1024
    elif block_size[-1].lower() == 'm':
        block_size = int(block_size[:-1]) * 1024 * 1024
    elif block_size[-1].lower() == 'k':
        block_size = int(block_size[:-1]) * 1024
    elif block_size[-1].lower() == 'b':
        block_size = int(block_size[:-1])
    else:
        block_size = int(block_size)

    sorter = ExternalSort(block_size)
    sorter.sort(filenames, os.path.join(temp_dir, outfilename))

    print('Outlier_DP threshold: {}'.format(get_quartile(os.path.join(temp_dir, 'dp.' + outfilename))))
    print('Outlier_GQ threshold: {}'.format(get_quartile(os.path.join(temp_dir, 'gq.' + outfilename))))
    os.remove(os.path.join(temp_dir, 'dp.' + outfilename))
    os.remove(os.path.join(temp_dir, 'gq.' + outfilename))

# for testing only
if __name__ == '__main__':
    filenames = sys.argv[1:]
    set_outlier(filenames, os.getcwd(), 'test.out', '500k')