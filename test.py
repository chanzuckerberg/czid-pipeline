import threading
import time

iostream = threading.Semaphore(2)


class AsyncHandler:
    def __init__(self):
        self.threads = []
        pass

    def launch(self, target, args):
        t = threading.Thread(target=target, args=args)
        self.threads.append(t)
        t.start()

    def wait_on_all(self):
        # write_to_log("Waiting on AsyncHandler threads...")
        for t in self.threads:
            t.join()
            print(t.result)
        # write_to_log("AsyncHandler threads finished.")

    def launch_aws_upload(self, src, dst):
        self.launch(self.aws_cp_work, (src, dst))

    def aws_cp_work(self, src, dst):
        with iostream:
            execute_command("aws s3 cp --quiet %s %s" % (src, dst))

    def launch_command(self, command):
        self.launch(execute_command, command)


def execute_command(cmd):
    print("Uploading file...: " + cmd)
    time.sleep(5)


async_handler = AsyncHandler()

async_handler.launch_aws_upload("input", "output")
async_handler.launch_aws_upload("output stuff lol", "out")
async_handler.launch_aws_upload("command 4", "out")
async_handler.launch_aws_upload("command 5", "out")
async_handler.wait_on_all()
print("everything is done!")
