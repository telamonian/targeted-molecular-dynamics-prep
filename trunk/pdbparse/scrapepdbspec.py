import urllib2
import BeautifulSoup
import re

Newlines = re.compile(unichr(160), flags=re.U)  #to deal with unicode whitespace weirdness

def GetPageText(url):
    # given a url, get page content
    data = urllib2.urlopen(url).read()
    # parse as html structured document
    bs = BeautifulSoup.BeautifulSoup(data, convertEntities=BeautifulSoup.BeautifulSoup.HTML_ENTITIES)
    # kill javascript content
    for s in bs.findAll('script'):
        s.replaceWith('')
    # find body and extract text
    txt = bs.find('body').getText('\n')
    #convert specific weird unicode whitespace to standard ascii whitespace
    return Newlines.sub(' ', txt)

def ScrapePDBSpec():
    urls = []
    html_template = 'http://www.wwpdb.org/documentation/format33/%s.html'
    for i in range (2,12):
        urls.append(html_template % ('sect'+str(i)))
    urls.append(html_template % 'remarks')
    print urls
    txt = ''.join([GetPageText(url) for url in urls])
    pdb_spec_scraped = open('pdb_spec_scraped.txt', 'w')
    pdb_spec_scraped.write(txt)

if __name__=="__main__":
    ScrapePDBSpec()