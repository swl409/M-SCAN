"""Utilities for parsing Human-GEM GPR rules used by M-SCAN."""

from __future__ import annotations

import math
import re


TOKEN_RE = re.compile(r"[A-Za-z0-9_]+|and|or|\(|\)")


class Parser:
    def __init__(self, tokens: list[str]):
        self.tokens = tokens
        self.pos = 0

    def peek(self) -> str | None:
        return self.tokens[self.pos] if self.pos < len(self.tokens) else None

    def consume(self, expected: str | None = None) -> str:
        token = self.peek()
        if token is None:
            raise ValueError("Unexpected end of GPR")
        if expected is not None and token != expected:
            raise ValueError(f"Expected {expected}, got {token}")
        self.pos += 1
        return token

    def parse(self):
        node = self.parse_or()
        if self.peek() is not None:
            raise ValueError(f"Unexpected token: {self.peek()}")
        return node

    def parse_or(self):
        nodes = [self.parse_and()]
        while self.peek() == "or":
            self.consume("or")
            nodes.append(self.parse_and())
        return flatten("or", nodes)

    def parse_and(self):
        nodes = [self.parse_factor()]
        while self.peek() == "and":
            self.consume("and")
            nodes.append(self.parse_factor())
        return flatten("and", nodes)

    def parse_factor(self):
        token = self.peek()
        if token == "(":
            self.consume("(")
            node = self.parse_or()
            self.consume(")")
            return node
        if token is None or token in {"and", "or", ")"}:
            raise ValueError(f"Unexpected token: {token}")
        self.consume()
        return ("gene", token)


def flatten(op: str, nodes: list):
    flat = []
    for node in nodes:
        if isinstance(node, tuple) and node[0] == op:
            flat.extend(node[1])
        else:
            flat.append(node)
    return flat[0] if len(flat) == 1 else (op, flat)


def parse_gpr(gpr: str):
    tokens = TOKEN_RE.findall(str(gpr))
    return Parser(tokens).parse()


def collect_genes(node) -> set[str]:
    if node is None:
        return set()
    if node[0] == "gene":
        return {node[1]}
    genes: set[str] = set()
    for child in node[1]:
        genes.update(collect_genes(child))
    return genes


def is_finite_value(value: float) -> bool:
    try:
        return not math.isnan(float(value))
    except (TypeError, ValueError):
        return False


def eval_gpr(node, gene_values: dict[str, float]) -> float:
    if node is None:
        return float("nan")
    kind = node[0]
    if kind == "gene":
        return float(gene_values.get(node[1], float("nan")))
    values = [eval_gpr(child, gene_values) for child in node[1]]
    values = [value for value in values if is_finite_value(value)]
    if not values:
        return float("nan")
    if kind == "or":
        return sum(values)
    if kind == "and":
        return min(values)
    raise ValueError(kind)
