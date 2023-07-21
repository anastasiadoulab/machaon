import sys
sys.path.append('..')
from src.grpc_core.jobreceiver import JobReceiver
from src.grpc_core import jobreceiver_pb2_grpc
from concurrent import futures
import grpc
import os

# Enviroment variables for debugging:
# os.environ['GRPC_VERBOSITY'] = 'debug'
# os.environ['GRPC_TRACE'] = 'all'

# The main class that handles gRPC communication
class gRPCServer:

    def __init__(self, port, root_disk, output_path, cores):
        self.port = port
        self.available_threads = 5
        self.server = grpc.server(futures.ThreadPoolExecutor(max_workers=self.available_threads))
        job_receiver = JobReceiver()
        job_receiver.setup_directories(root_disk, output_path, cores)
        jobreceiver_pb2_grpc.add_JobReceiverServicer_to_server(job_receiver, self.server)
        certs_path =  '<CERTS_PATH>'
        with open(os.path.join(certs_path, 'node1.key'), 'rb') as key_file:
            private_key = key_file.read()
        with open(os.path.join(certs_path, 'node1.cert'), 'rb') as cert_file:
            cert_chain = cert_file.read()
        with open(os.path.join(certs_path, 'machaonlocalca.cert'), 'rb') as root_cert_file:
            root_cert = root_cert_file.read()
        credentials = grpc.ssl_server_credentials([(private_key, cert_chain)], root_cert,
        require_client_auth=True)
        self.server.add_secure_port('[::]:' + self.port, credentials) 
        
    def start(self):
        self.server.start()
        print("Server started, listening on " + self.port)
        self.server.wait_for_termination()

if __name__ == '__main__':
    # grpc_server = gRPCServer('8080', '/opt/storage/machaon_data',
    #                          '/opt/storage/output', 8)
    grpc_server = gRPCServer('55555', '<DATA_DIR>',
                             '<OUTPUT_DIR>', 30)
    grpc_server.start()
    
   
