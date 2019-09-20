# cramore

`cramore` is a collection of C++ tools to manipulate SAM/BAM/CRAM and
BCF/VCF files in various contexts of sequence analysis.

### Installing cramore

Before installing `cramore`, you need to install
[htslib](https://github.com/samtools/htslib) in the same directory you
want to install demuxlet (i.e. demuxlet and htslib should be
siblings). You also need [cmake](https://cmake.org/) installed in your system.

After installing htslib, you can clone the current snapshot of this repository to install as well

<pre>
$ mkdir build

$ cd build

$ cmake ..
</pre>

In case any required libraries is missing, you may specify customized installing path by replacing "cmake .." with:

<pre>
For libhts:
  - $ cmake -DHTS_INCLUDE_DIRS=/hts_absolute_path/include/  -DHTS_LIBRARIES=/hts_absolute_path/lib/libhts.a ..

For bzip2:
  - $ cmake -DBZIP2_INCLUDE_DIRS=/bzip2_absolute_path/include/ -DBZIP2_LIBRARIES=/bzip2_absolute_path/lib/libbz2.a ..

For lzma:
  - $ cmake -DLZMA_INCLUDE_DIRS=/lzma_absolute_path/include/ -DLZMA_LIBRARIES=/lzma_absolute_path/lib/liblzma.a ..
</pre>

Finally, to build the binary, run

<pre>
$ make
</pre>

### List of tools contained in `cramore`

`cramore` contains many in-house C++ tools that are currently under
the hood development phase. To list the available commands of tools, type:

<pre>
cramore --help
</pre>

To see the usage of individual commands, type:

<pre>
cramore [command] --help
</pre>
