#!/usr/bin/python
import glob
import os
import time

def get_files(dirs):
  file_pairs=[]
  for root_dir in dirs:
    for root, directories, files in os.walk(root_dir):
      for filename in files:
        if (filename[-2:] == 'pp'):
          file_pairs.append([root,filename])
  return file_pairs


def file_stats(file_pairs):
  loc = 0
  nfiles = 0
  nsuites = 0
  ntests = 0
  for path, filename in file_pairs:
    loc += int(os.popen('wc -l '+path+'/'+filename).read().split()[0])
    nfiles+=1
    if (filename[:4] == 'Test'):
      nsuites+=1
      ntests+= int(os.popen('egrep -c -i "void\s+Test" '+path+'/'+filename).read().split()[0])
  return [nfiles, loc, nsuites, ntests]  

def print_stats():
  rev1_epoch=1113991642
  svn_info = os.popen('svn info').read().split('\n')
  rev_line=svn_info[4]
  revision = rev_line.split()[1]
  date_line = svn_info[9].split()
  date_time=date_line[3]+' '+date_line[4]
  pattern = '%Y-%m-%d %H:%M:%S'
  epoch = int(time.mktime(time.strptime(date_time, pattern))) - rev1_epoch
  
  
  source_dirs=glob.glob('*/src')
  test_dirs=glob.glob('*/test')+glob.glob('*/tests')
  
  test_files=get_files(test_dirs)
  source_files=get_files(source_dirs)

  test_stats=file_stats(test_files)
  source_stats=file_stats(source_files)

  print revision,'\t',epoch,'\t',source_stats[0],'\t',source_stats[1],'\t',test_stats[0],'\t',test_stats[1],'\t',test_stats[2],'\t',test_stats[3]


last_revision=3578
step=10

print '#rev\ttime\tsrc_files\tsrc_loc\ttest_files\ttests_loc\ttest_suites\ttests'
for rev in range(1,last_revision):
  if (rev%step == 0):
    os.system('svn up -r '+str(rev)+' > /dev/null')
    print_stats()
