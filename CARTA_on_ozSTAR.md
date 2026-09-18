# Running CARTA Remotely on OzStar

CARTA (the radio-astronomy image viewer) runs as a backend server on an OzStar login node, but that node has no GUI browser installed — so you can't view CARTA there directly. The fix: start the CARTA backend on the cluster, then open an SSH tunnel from your own laptop so you can view it in a browser on your own machine. Use this whenever you need to inspect FITS images or cubes on OzStar with CARTA.

## Step 1: Launch CARTA on the cluster

SSH into OzStar as usual, then:

```bash
tmux new -s carta
cd carta_directory
APPIMAGE_EXTRACT_AND_RUN=1 ./carta-x86_64.AppImage --no_browser
```

`tmux` keeps CARTA running if the SSH session drops — detach with `Ctrl+b` then `d`. `--no_browser` skips CARTA's harmless-but-noisy attempt to auto-launch a GUI browser that doesn't exist on the login node.

CARTA will print a line like:

```
CARTA is accessible at http://192.168.44.13:3002/?token=8001340d-c9fc-48a1-96b6-d834def947de
```

Note the port (usually 3002) and the token — you need both in Step 3.

## Step 2: Open an SSH tunnel from your own computer

In a **new terminal on your own laptop** (not the cluster), tunnel a local port to the internal cluster address CARTA printed:

```bash
ssh -N -L 3002:192.168.44.13:3002 username@ozstar.hpc.swin.edu.au
```

Replace the hostname with whatever you normally use to SSH into OzStar, and match both ports to whatever CARTA actually printed if it differs from 3002. `-N` just holds the tunnel open with no remote command — leave this terminal running for as long as you want to use CARTA.

## Step 3: Open CARTA in your local browser

With the tunnel open, go to:

```
http://localhost:3002/?token=8001340d-c9fc-48a1-96b6-d834def947de
```

using `localhost` in place of the internal cluster IP, but keeping the same port and token CARTA printed in Step 1. The CARTA UI should load normally in your own browser.

## Troubleshooting

- **"Loading user-defined snippets failed!" on first launch** — harmless. CARTA looks for a saved scripting-console snippets file in its config directory and errors instead of silently skipping it when none exists yet (typical on a first run on a new account/machine). Click OK and continue; CARTA works normally and creates the file once you save a snippet.
- **"CARTA Usage Data" prompt** — a one-time, unrelated opt-in for anonymous telemetry (session duration, image count/size). Either choice is fine; it doesn't affect functionality.
- **SSH tunnel refuses to connect** — you're probably not using the exact hostname/login node you normally SSH in with (e.g. `farnarkle1` vs `farnarkle2` vs a load-balanced name). Match it exactly.
- **Port or token mismatch** — always copy the exact port and token CARTA printed for *that specific run*; both can change between runs.
- **Session dies when you close your terminal** — always launch CARTA inside `tmux` (or `screen`) on the cluster side so it survives a dropped SSH connection.
