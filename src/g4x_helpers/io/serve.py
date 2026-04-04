#!/usr/bin/env python3

import argparse
import os
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer


class CORSRequestHandler(SimpleHTTPRequestHandler):
    def end_headers(self):
        self.send_header('Access-Control-Allow-Origin', '*')
        self.send_header('Access-Control-Allow-Methods', 'GET, HEAD, OPTIONS')
        self.send_header('Access-Control-Allow-Headers', '*')
        super().end_headers()

    def do_OPTIONS(self):
        self.send_response(204)
        self.end_headers()


def serve(directory: str, port: int = 8000):
    os.chdir(directory)

    httpd = ThreadingHTTPServer(('0.0.0.0', port), CORSRequestHandler)

    host, port = httpd.server_address
    print(f'Serving directory: {os.getcwd()}')
    print(f'URL: http://{host}:{port}')

    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        print('\nShutting down server...')
        httpd.server_close()


def main():
    parser = argparse.ArgumentParser(description='Simple HTTP server with CORS support')
    parser.add_argument('-d', '--directory', default='.', help='Directory to serve (default: current directory)')
    parser.add_argument('-p', '--port', type=int, default=8000, help='Port to serve on (default: 8000)')

    args = parser.parse_args()

    serve(args.directory, args.port)


if __name__ == '__main__':
    main()
