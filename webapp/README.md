# TTQ webapp

This folder contains the webapp. The webapp consists of a Jupyter notebook which is served via [Voilà](https://voila.readthedocs.io/).

To install, run `pip install -e .[web]`. If everything is working, `./start_server` should start the Voilà server.

## Detailed Linux VM setup instructions

This assumes you have Anaconda or [Miniconda](https://docs.conda.io/en/latest/miniconda.html) installed.

1. Clone this repository to the VM
2. Create a new screen session: `screen -S ttq_app`
3. Create a new Conda environment: `conda create -n ttq_app python=3.8`
4. Activate: `conda activate ttq_app`
5. Install main dependencies with `pip install -e .[web]` in the main folder
6. Edit `server_name` in `nginx_ttq` to match your local server
7. Symlink `nginx_ttq`  to `/etc/nginx/sites-enabled`
8. Restart NGINX: `sudo service nginx reload`
9. Start the app: `./start_server` in this folder
10. If/when desired, disconnect from the screen session: `ctrl+a` then `d`