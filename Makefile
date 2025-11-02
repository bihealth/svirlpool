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
	@echo "  lock           Build uv.lock file"

.PHONY: check
check:
	uv run hatch run quality:check
	uv run hatch run quality:typecheck

.PHONY: fix
fix:
	uv run hatch run quality:format

.PHONY: test
test:
	uv run hatch run tests:run

.PHONY: test-snapshot
test-snapshot:
	uv run hatch run tests:run-snapshot

.PHONY: lock
lock:
	uv lock --upgrade
