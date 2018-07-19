# cramore

`cramore` is a collection of C++ tools to manipulate SAM/BAM/CRAM and
BCF/VCF files in various contexts of sequence analysis.

### Installing `cramore`

Before installing `cramore`, you need to install [htslib](https://github.com/samtools/htslib) in the same directory you want to install demuxlet (i.e. demuxlet and htslib should be siblings).

After installing htslib, you can clone the current snapshot of this repository to install as well

<pre>
$ git clone https://github.com/hyunminkang/cramore.git
$ cd cramore
$ autoreconf -vfi
$ ./configure (with additional options such as --prefix)
$ make
$ make install (may require root privilege depending on the prefix)
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
