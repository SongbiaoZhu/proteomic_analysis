# -*- coding: utf-8 -*-
# Author: fcj
# Time: 2019-05-09
# Description: python 多线程-普通多线程-生产者消费者模型

import os
import time
import requests
import threading
from bs4 import BeautifulSoup
from queue import Queue
from threading import Thread

# extract the protein id lines  
org = 'Rat'
orgFull = 'RAT'

publicDir = 'E:/publicData/' 
outDir = publicDir + 'uniprotFunctions{}/'.format(org)

def uniFuncSpider(url):
    headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/59.0.3071.104 Safari/537.36'
    }
    r = requests.get(url=url, headers=headers)
    r.raise_for_status()
    r.encoding = r.apparent_encoding
    html = r.text
    soup = BeautifulSoup(html, "lxml")
    try:
        x = soup.select('#function > div:nth-child(2) > span')[0]
        func = x.get_text().split('\r\n')[0]                   
    except:
        func = 'Not found function summary.'
    fname = url.split('/')[-1] + '.txt'
    fw = open(os.path.join(outDir, fname), 'w')
    fw.write(func)
    fw.close()
    return None

def producer(in_q):  # 生产者
    ready_list = []
    base_url = 'https://www.uniprot.org/uniprot/'
    # 构造所有url
    f = open(os.path.join(publicDir, 'Accession{}.txt'.format(org)), 'r')
    accession = f.readlines()
    f.close()    
    while in_q.full() is False:
        for i in accession[7971:]:
            url = base_url + i.replace("\n", "")
            if url not in ready_list:
                ready_list.append(url)
                in_q.put(url)
            else:
                continue
'''
Check which errored Accessions, craw again
crawedAccession = []
for root, dirs, files in os.walk(outDir):
    for i in files:
        crawedAccession.append(i.split('.')[0])
f = open(os.path.join(publicDir, 'Accession{}.txt'.format(org)), 'r')
allAccession = f.readlines()
f.close()   
newAccession = [i.strip() for i in allAccession[1:]]
list0 = list(set(newAccession)-set(crawedAccession)) 

def producer(in_q):  # 生产者
    ready_list = []
    base_url = 'https://www.uniprot.org/uniprot/'  
    while in_q.full() is False:
        for i in list0:
            url = base_url + i.replace("\n", "")
            if url not in ready_list:
                ready_list.append(url)
                in_q.put(url)
            else:
                continue
'''

def consumer(in_q, out_q):  # 消费者
    while True:
        url =in_q.get()
        uniFuncSpider(url)
        out_q.put(str(threading.current_thread().getName()))
        in_q.task_done()  # 通知生产者，队列已消化完

if __name__ == '__main__':
    start = time.time()
    if not os.path.exists(outDir):
        os.mkdir(outDir)    
    queue = Queue(maxsize=40)  # 设置队列最大空间为10
    result_queue = Queue()
    print('queue 开始大小 %d' % queue.qsize())

    producer_thread = Thread(target=producer, args=(queue,))
    producer_thread.daemon = True
    producer_thread.start()

    for index in range(40):
        consumer_thread = Thread(target=consumer, args=(queue, result_queue, ))
        consumer_thread.daemon = True
        consumer_thread.start()

    queue.join()
    end = time.time()
    print('总耗时：%s' % (end - start))
    print('queue 结束大小 %d' % queue.qsize())
    print('result_queue 结束大小 %d' % result_queue.qsize())
