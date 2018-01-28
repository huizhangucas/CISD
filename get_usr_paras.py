# -*- coding:utf8 -*-
import os
import datetime
import sys
from optparse import OptionParser


def get_user_paras():
    try:
        opt = OptionParser()
        opt.add_option('--host_ip',
                       dest='host_ip',
                       type=str,
                       help='the ip of the check host')
        opt.add_option('--run',
                       action="store_true",
                       dest="is_run",
                       default=False,
                       help="run the scripts")
        opt.add_option('--view',
                       action="store_false",
                       dest="is_run",
                       default=False,
                       help="only view but not run the scripts")
        opt.add_option('--show_type',
                       dest="show_type",
                       type=int,
                       default=0,
                       help="0 or 1, 0 only show the simple data, 1 show the full data")
        (options, args) = opt.parse_args()
        is_valid_paras = True
        error_messages = []
        host_ip = options.host_ip
        is_run = options.is_run
        show_type = options.show_type
        if not host_ip:
            error_messages.append("host_ip must be set;")
            is_valid_paras = False
        if show_type not in [0, 1]:
            error_messages.append("show_type only can be 0 or 1;")
            is_valid_paras = False

        if is_valid_paras:
            user_paras = {"host_ip": host_ip, "is_run": is_run, "show_type": show_type}
            return user_paras
        else:
            for error_message in error_messages:
                print(error_message)
               opt.print_help()
            return None
    except Exception as ex:
        print("exception :{0}".format(str(ex)))
        return None


def main():
    user_paras = get_user_paras()
    if user_paras is None:
        sys.exit(0)
    info = "host_ip:{0}, is_run:{1}, show_type:{2}"
    info = info.format(user_paras["host_ip"],
                       user_paras["is_run"],
                       user_paras["show_type"])
    print(info)


if __name__ == '__main__':
    main()