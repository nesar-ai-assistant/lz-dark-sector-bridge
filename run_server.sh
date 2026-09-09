#!/bin/bash
# Wrapper script to run the LZ Dark Sector Bridge MCP server
# Ensures correct working directory and Python environment
cd "$(dirname "$0")"
exec python -m mcp_server "$@"
