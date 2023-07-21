import grpc
import jobreceiver_pb2
import jobreceiver_pb2_grpc

def run():
    with grpc.insecure_channel('localhost:50505') as channel:
        stub = jobreceiver_pb2_grpc.JobReceiverStub(channel)
        request = jobreceiver_pb2.ResultRequest(request_id=9, hash='9955620089909220348')
        responses = stub.DownloadResult(request)
        for response in responses:
            if response.WhichOneof('job_data') == 'file_info':
                print(f"Received job details: {response.file_info}")
            elif response.WhichOneof('job_data') == 'chunk_data':
                # Write the chunk data to file or process it in some way
                print(f"Received chunk data of size {len(response.chunk_data)}")


if __name__ == '__main__':
    run()