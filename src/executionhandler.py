import tqdm
from pebble import ProcessPool
import time
import sys
import os

class ExecutionHandler:

    def __init__(self, cores, max_seconds):
        self.progress_bar = None
        self.max_wait_time = max_seconds
        self.cores = cores
        self.chunk_size = 10000

    def parallelize(self, func, entries):
        # Run processes in parallel with a timeout
        start = time.time()
        count = 0
        results = []
        chunks = int(len(entries) / self.chunk_size)
        if (len(entries) % self.chunk_size != 0):
            chunks += 1
        with tqdm.tqdm(total=len(entries), position=0, leave=True, file=sys.stdout) as progress_bar:
            for chunk in range(0, chunks):
                slice_start = self.chunk_size * chunk
                slice_end = self.chunk_size * (chunk+1)
                input_chunk = entries[slice_start:slice_end]
                with ProcessPool(max_workers=self.cores) as pool:
                    future = pool.map(func, input_chunk, timeout=self.max_wait_time)
                    iterator = future.result()
                    while True:
                        count += 1
                        try:
                            result = next(iterator)
                            results.append(result)
                            progress_bar.update(1)
                        except StopIteration:
                            break
                        except TimeoutError as error:
                            print("job %d took longer than %d seconds" % (count, error.args[1]))
                        except Exception as error:
                            print("job %d  raised %s" % (count, error))
                            print(repr(error))
        end = time.time()
        total = len(results)
        final_results = [result for result in results if result is not False]
        if(total != len(final_results)):
            fail_log = [repr(entries[resultIndex]) for resultIndex, result in enumerate(results) if result is False]
            with open(os.path.join('logs', ''.join([time.strftime('%Y-%m-%d', time.localtime()), '.log'])), 'a', encoding='utf-8') as fail_log_file:
                fail_log_file.write('\n'.join([repr(func)] + fail_log + ['\n']))
        sys.stdout.flush()
        print(''.join(['Execution time (seconds): ', repr(round(end - start, 3)), ' | fails: ', repr(total - len(final_results))]))
        return final_results
