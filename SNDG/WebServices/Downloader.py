from SNDG.WebServices import download_file
class FileDownload():

    def __init__(self,url,filename=False):
        self.url = url
        self.filename = filename if filename else url.split("/")[-1]
        self.format = "flat"
        self._uncompress = True

    def download(self,outdir):
        download_file

    def override(self):
        pass

    def check_old(self,days=7):
        pass

    def min_length(self,min_len:int):
        self.min_len = min_len
        return self
    def unzip(self):
        self.format = "zip"
        return self
    def gunzip(self):
        self.format = "gz"
        return self
    def donot_uncompress(self):
        self._uncompress = False
        return self
    def md5(self,md5):
        self._md5 = md5
        return self
    def download_md5(self,md5_url):
        self._download_md5 = md5_url
        return self

class FolderDownload():
    pass

class MixDownload():
    pass


class Downloader():

    def __init__(self):
        pass


