from distutils.core import setup

__version__ = '0.0.1'

setup(name='rfipip',
      version=__version__,
      description='',
      long_description='',
      license='GPL',
      author='Monika Obrocka',
      author_email='mobrocka at ska.ac.za',
      classifiers=[
          'Development Status :: 3 - Alpha',
          'Intended Audience :: Developers',
          'Operating System :: OS Independent',
          'License :: OSI Approved :: GNU General Public License (GPL)',
          'Topic :: Scientific/Engineering :: Astronomy',
          'Topic :: Software Development :: Libraries :: Python Modules',
      ],
      requires=['h5py', 'sqlite3', 'matplotlib', 'numpy', 'csv', 'pylab', 'sigpyproc'],
      provides=['rfipip'],
      package_dir={'rfipipeline': 'rfipip'},
      packages=['rfipip'],
      url='https://github.com/pinsleepe/rfipip.git',  # use the URL to the github repo
      download_url='https://github.com/pinsleepe/rfipip/archive/0.0.1.tar.gz',
      keywords=['rfi', 'astronomy', 'ml'],  # arbitrary keywords
      )
