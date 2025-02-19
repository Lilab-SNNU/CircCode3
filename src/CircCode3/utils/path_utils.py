# -*- coding = utf-8 -*-
import os
import time
import shutil

from ..utils.logs import logger


class PathManager:
    def __init__(self, parsent = {}):
        self.paths: dict[str, list[str]] = parsent
    
    def __getitem__(self, item):
        return self.paths[item]

    def __setitem__(self, key, value):
        self.paths[key] = value

    def add_batch_paths(self, paths: dict[str, str]):
        self.paths = {**self.paths, **paths}

    def add_path(self, keyword, path):
        self.paths[keyword] = path

    def get_path(self, keyword):
        return self.paths.get(keyword, None)

    def remove_path(self, keyword, path):
        if keyword in self.paths:
            self.paths[keyword].remove(path)
            if not self.paths[keyword]:
                self.paths.pop(keyword)

def create_folder(path):
    if not os.path.exists(path):
        os.mkdir(path)


def del_all_file(remove_path):
    if not os.path.exists(remove_path):
        return
    if os.path.isfile(remove_path):
        try:
            os.remove(remove_path)
        except BaseException as e:
            logger.error(f"Error deleting temporary file: \n{e}")
    elif os.path.isdir(remove_path):
        file_lis = os.listdir(remove_path)
        for file_name in file_lis:
            tf = os.path.join(remove_path, file_name)
            del_all_file(tf)
        if not os.listdir(remove_path):
            os.rmdir(remove_path)


def copyfile(source_file, arm_path):
    fpath, fname = os.path.split(source_file)
    shutil.copy(source_file, arm_path + "/" + fname)
    return arm_path + "/" + fname


def project_name(outpath):
    # 生成项目名
    project = outpath + '/project'
    now = time.strftime('%y-%m-%d %H:%M:%S', time.localtime())
    projectname = project + now.split()[0]  # 项目名“project” + 日期
    num = 1
    while os.path.exists(projectname + f'.{num}'):
        num += 1
    projectname = projectname + f'.{num}'
    return projectname

