#!/usr/bin/env python3

import json
import sys

print("loading stdin")
data = json.load(sys.stdin)
print(f"loading {sys.argv[1]}")
reference = json.load(open(sys.argv[1]))
print("checking...")
if data != reference:
    print("Result does not match reference.")
    #print(f"reference={reference}")
    #print(f"data={data}")
    sys.exit(1)
else:
    print("Correct!")
    sys.exit(0)
    
