.PHONY: default
default: help

.PHONY: help
help:
	@echo "Usage: make [target]"
	@echo ""
	@echo "Targets:"
	@echo "  help       	Show this help message"
	@echo "  check      	Check the project"
	@echo "  fix        	Fix the project"
	@echo "  test       	Run the tests"
	@echo "  test-snapshot	Run the tests and rebuild snapshots"
	@echo "  lock           Update pixi.lock file"

.PHONY: check
check:
	pixi run -e dev check
	pixi run -e dev typecheck

.PHONY: fix
fix:
	pixi run -e dev format

.PHONY: test
test:
	pixi run -e dev test

.PHONY: test-snapshot
test-snapshot:
	pixi run -e dev test-snapshot

.PHONY: lock
lock:
	pixi update
