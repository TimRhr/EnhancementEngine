# Web App Guide

This short guide shows how to start the optional Flask web interface.

## Installation

Install Flask if it is not already available:

```bash
pip install Flask
```

## Running the App

A small example application is provided in the `examples` folder. Start it with:

```bash
python examples/webapp.py
```

Once the server is running, open `http://127.0.0.1:5000/` in your browser to
access the form-based interface. You can also send a POST request to `/analyze`
containing `gene` and `variant` fields to receive a JSON analysis report.

