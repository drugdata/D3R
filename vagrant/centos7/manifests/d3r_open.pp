class d3r_open
{
  Package { ensure => 'installed' }

  $python_deps    = [ 'epel-release', 'python-pip', 'python-psutil', 'python-virtualenv', 'python-tox', 'pylint', 'python-coverage' ]
  $perl_deps      = [ 'perl-Archive-Tar', 'perl-List-MoreUtils' ]
  $other_packages = [ 'libXft', 'openbabel', 'xorg-x11-xauth' ]
  $pymol_deps     = [ 'subversion', 'gcc', 'gcc-c++', 'kernel-devel', 'python-devel', 'tkinter', 'python-pmw', 'glew-devel', 'freeglut-devel', 'libpng-devel', 'freetype-devel', 'libxml2-devel']
  $mesa_packages  = [ 'mesa-libGL-devel','mesa-libEGL-devel','mesa-libGLES-devel' ]
  $pip_packages   = [ 'argparse','psutil','biopython','xlsxwriter','ftpretty','wheel','flake8','lockfile','easywebdav' ]
  package { $python_deps: }
  package { $perl_deps: }
  package { $other_packages: }
  package { $mesa_packages: }
  package { $pymol_deps: }

  #package { "PyMOL": 
  #  ensure => installed 
  #}

  exec { 'install_pymol':
    command => '/bin/cd ~; /bin/wget https://heanet.dl.sourceforge.net/project/pymol/pymol/1.8/pymol-v1.8.4.0.tar.bz2;
                /bin/tar -jxvf  pymol-v1.8.4.0.tar.bz2;
                cd pymol; 
                python2.7 setup.py build install --home=/opt/pymol --install-scripts=/opt/pymol --install-lib=/opt/pymol/modules;
                /bin/ln -s /opt/pymol/pymol /usr/bin/pymol',

    creates => '/opt/pymol/pymol'
  }
  yumrepo { 'giallu-rdkit':
    baseurl             => 'https://copr-be.cloud.fedoraproject.org/results/giallu/rdkit/epel-$releasever-$basearch/',
    descr               => 'Copr repo for rdkit owned by giallu',
    enabled             => 1,
    gpgcheck            => 1,
    gpgkey              => 'https://copr-be.cloud.fedoraproject.org/results/giallu/rdkit/pubkey.gpg',
    repo_gpgcheck       => 0,
    skip_if_unavailable => "True"
  }

  exec { 'install_conda':
    command => '/bin/cd ~; /bin/wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh;
                /bin/chmod a+x Miniconda2-latest-Linux-x86_64.sh;
                /bin/sh Miniconda2-latest-Linux-x86_64.sh -b -p /opt/miniconda2;
                export PATH="/opt/miniconda2/bin:$PATH";
                conda update --yes conda;
                conda install -y -c rdkit rdkit=2016.03.3',
    creates => '/opt/miniconda2'
  }

  package { 'python-rdkit':
    require => Yumrepo['giallu-rdkit']
  }
  package { $pip_packages:
    ensure   => 'installed',
    provider => 'pip',
    require  => Package['python-pip']
  }
  #blast
  exec { 'install_blast':
    command => '/bin/yum install -y ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.3.0/ncbi-blast-2.3.0+-1.x86_64.rpm',
    creates => '/usr/bin/blastp'
  }
  #autodock
  exec { 'install_autodock':
    command => '/bin/wget http://vina.scripps.edu/download/autodock_vina_1_1_2_linux_x86.tgz;
                /bin/tar -zxf autodock_vina_1_1_2_linux_x86.tgz -C /opt;
                /bin/ln -s /opt/autodock_vina_1_1_2_linux_x86/bin/vina /usr/local/bin/vina; 
                /bin/ln -s /opt/autodock_vina_1_1_2_linux_x86/bin/vina_split /usr/local/bin/vina_split',
    creates => '/opt/autodock_vina_1_1_2_linux_x86/bin/vina'
  }
  exec { 'install_mgltools':
    command => '/bin/wget http://mgltools.scripps.edu/downloads/downloads/tars/releases/REL1.5.6/mgltools_x86_64Linux2_1.5.6.tar.gz;
                /bin/tar -zxf mgltools_x86_64Linux2_1.5.6.tar.gz;
                cd mgltools_x86_64Linux2_1.5.6; 
                /bin/sh install.sh -d /usr/local/mgltools;
                export MGL_ROOT=/usr/local/mgltools >> ~/.bash_profile',
    creates => '/usr/local/mgltools'
  }

  #manual Install UCSF chimera
  #manual INSTALL Schrodinger_Suites_2016-2_Linux-x86_64

}

# Run the d3r_open class
class { 'd3r_open': }
