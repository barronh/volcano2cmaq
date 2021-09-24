import argparse

if __name__ == '__main__':
    from . import runcfg
    parser = argparse.ArgumentParser()
    parser.add_argument('config', nargs='*', default=['run.cfg'])

    args = parser.parse_args()
    runcfg(args.config)
