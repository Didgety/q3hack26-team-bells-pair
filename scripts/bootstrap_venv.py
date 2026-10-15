"""Create this project's isolated Python environment; optionally install inputs."""
import argparse
from pathlib import Path
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--install', action='store_true',
                        help='Install requirements.in using the .venv interpreter.')
    parser.add_argument('--locked', action='store_true',
                        help='Install the tested requirements.lock.txt instead of resolving inputs.')
    parser.add_argument('--python', default=sys.executable,
                        help='CPython executable used to create .venv; defaults to this interpreter.')
    args = parser.parse_args()
    root = Path(__file__).resolve().parent.parent
    env = root / '.venv'
    env_python = env / ('Scripts/python.exe' if sys.platform == 'win32' else 'bin/python')
    if not env.exists():
        subprocess.run([args.python, '-m', 'venv', str(env)], check=True)
    if not env_python.is_file():
        raise SystemExit(f'Existing .venv has no expected interpreter: {env_python}')
    subprocess.run([str(env_python), '-c',
                    "import sys; assert sys.prefix != sys.base_prefix; "
                    "print('Virtualenv interpreter:', sys.executable); "
                    "print('Python:', sys.version.split()[0])"], check=True)
    if args.install:
        requirements = root / ('requirements.lock.txt' if args.locked else 'requirements.in')
        if not requirements.is_file():
            raise SystemExit(f'Missing requirements file: {requirements}')
        subprocess.run([str(env_python), '-m', 'pip', 'install',
                        '-r', str(requirements)], check=True)
        subprocess.run([str(env_python), '-m', 'pip', 'check'], check=True)
    print('Use this interpreter for scripts, tests, pip, and notebook kernels:', env_python)


if __name__ == '__main__':
    main()
