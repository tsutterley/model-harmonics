jpl_ecco_webdav.py
==================

- Retrieves and prints a user's ECCO Drive WebDAV credentials
- WebDAV credentials can be used in the [ECCO Drive sync program](https://github.com/tsutterley/model-harmonics/blob/main/ECCO/jpl_ecco_sync.py)

#### Calling Sequence
```bash
python jpl_ecco_webdav.py --directory <path_to_directory>
```
[Source code](https://github.com/tsutterley/model-harmonics/blob/main/ECCO/jpl_ecco_webdav.py)

#### Command Line Options
- `-U X`, `--user X`: Username for NASA Earthdata Login
- `-N X`, `--netrc X`: Path to .netrc file for authentication
- `-A`, `--append`: Append .netrc file instead of printing
