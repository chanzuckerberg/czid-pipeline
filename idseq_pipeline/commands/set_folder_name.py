"""Pre-processing commands."""
from .base import Base


class SetFolderName(Base):
    """Perform pre-processing functions"""

    def run(self):
        self.set_index_date(['/blast/db/FASTA/nr.gz', '/blast/db/FASTA/nt.gz'])
        pass

    def set_index_date(self, input_files):
        from .common import *

        # install_ncbitool('.')
        ncbitool_path = './ncbitool'
        maxe = None
        for f in input_files:
            command = "%s file %s" % (ncbitool_path, f)
            output = execute_command_with_output(command).split("File Info: ")[1]
            date = json.loads(output)["ModTime"]
            print date
            date = datetime.datetime.strptime(date, "%Y-%m-%dT%H:%M:%S")
            if maxe is None or date > maxe:
                maxe = date
            print date
        print "MAX:"
        print maxe

        source_time = datetime.datetime.strftime(maxe, "%Y-%m-%d-utc-")
        source_utime = str(int(time.mktime(maxe.timetuple())))
        source_str = source_time + source_utime + "-unixtime"
        print(source_str)

        now = datetime.datetime.utcnow()
        utime = str(int(time.time()))
        now_str = datetime.datetime.strftime(now, "%Y-%m-%d-utc-") + utime + "-unixtime"
        print(now_str)
        print("Final date:")
        print(source_str + "__" + now_str)
