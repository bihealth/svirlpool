# Svirlpool: structural variant detection from long read sequencing by local assembly

## Installing & Running

### Prerequisites

### Instructions

## Developer Notes

### Prerequisites

Install [`pixi`](https://pixi.sh/latest/):

```
curl -fsSL https://pixi.sh/install.sh | bash
```

### Instructions

Open VS Code with pix-installed `dev` environment

```
pixi run -e dev code .
```

Clone repository:

```
git clone git@github.com:bihealth/svirlpool.git
cd svirlpool
```

Run code formatting

```
make fix
```

Run formatting, lints, other static checks:

```
make check
```

Run tests:

```
make check
```

Or all in one:

```
make fix check test
```
