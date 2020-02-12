// Copyright © 2018-2020 Wei Shen <shenwei356@gmail.com>
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

package unikmer

import (
	"errors"
	"fmt"
	"strconv"
	"strings"

	"github.com/shenwei356/breader"
)

// Taxonomy holds relationship of taxon in a taxonomy.
type Taxonomy struct {
	file     string
	rootNode uint32

	Nodes map[uint32]uint32 // parent -> child

	cacheLCA bool
	lcaCache map[uint64]uint32 // cache of lca

	maxTaxid uint32
}

// ErrIllegalColumnIndex means column index is 0 or negative.
var ErrIllegalColumnIndex = errors.New("unikmer: illegal column index, positive integer needed")

// NewTaxonomyFromNCBI parses Taxonomy from nodes.dmp
// from ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz .
func NewTaxonomyFromNCBI(file string) (*Taxonomy, error) {
	return NewTaxonomy(file, 1, 3)
}

// NewTaxonomy loads nodes from file,
// for example,
func NewTaxonomy(file string, childColumn int, parentColumn int) (*Taxonomy, error) {
	if childColumn < 1 || parentColumn < 1 {
		return nil, ErrIllegalColumnIndex
	}
	minColumns := childColumn
	if parentColumn > minColumns {
		minColumns = parentColumn
	}

	// taxon represents a taxonomic node
	type taxon struct {
		Taxid  uint32
		Parent uint32
	}

	parseFunc := func(line string) (interface{}, bool, error) {
		items := strings.Split(line, "\t")
		if len(items) < minColumns {
			return nil, false, nil
		}
		child, e := strconv.Atoi(items[childColumn-1])
		if e != nil {
			return nil, false, e
		}
		parent, e := strconv.Atoi(items[parentColumn-1])
		if e != nil {
			return nil, false, e
		}
		return taxon{Taxid: uint32(child), Parent: uint32(parent)}, true, nil
	}

	reader, err := breader.NewBufferedReader(file, 8, 100, parseFunc)
	if err != nil {
		return nil, fmt.Errorf("unikmer: %s", err)
	}

	nodes := make(map[uint32]uint32, 1024)
	var root uint32

	var tax taxon
	var data interface{}
	var maxTaxid uint32
	for chunk := range reader.Ch {
		if chunk.Err != nil {
			return nil, fmt.Errorf("unikmer: %s", err)
		}
		for _, data = range chunk.Data {
			tax = data.(taxon)

			nodes[tax.Taxid] = tax.Parent

			if tax.Taxid == tax.Parent {
				root = tax.Taxid
			}
			if tax.Taxid > maxTaxid {
				maxTaxid = tax.Taxid
			}
		}
	}

	return &Taxonomy{file: file, Nodes: nodes, rootNode: root, maxTaxid: maxTaxid}, nil
}

// MaxTaxid returns maximum taxid
func (t *Taxonomy) MaxTaxid() uint32 {
	return t.maxTaxid
}

// CacheLCA tells to cache every LCA query result
func (t *Taxonomy) CacheLCA() {
	t.cacheLCA = true
	if t.lcaCache == nil {
		t.lcaCache = make(map[uint64]uint32, 1024)
	}
}

// LCA returns the Lowest Common Ancestor of two nodes
func (t *Taxonomy) LCA(a uint32, b uint32) uint32 {
	if a == 0 {
		return b
	}
	if b == 0 {
		return a
	}
	if a == b {
		return a
	}

	// check cache
	var ok bool

	var query uint64
	if t.cacheLCA {
		query = pack2uint32(a, b)
		var c uint32
		if c, ok = t.lcaCache[query]; ok {
			return c
		}
	}

	lineA := make([]uint32, 0, 16)
	mA := make(map[uint32]struct{}, 16)

	var child, parent uint32

	child = a
	for {
		parent, ok = t.Nodes[child]
		if !ok {
			t.lcaCache[query] = b
			return 0
		}
		if parent == child { // root
			lineA = append(lineA, parent)
			mA[parent] = struct{}{}
			break
		}
		if parent == b { // b is ancestor of a
			if t.cacheLCA {
				t.lcaCache[query] = b
			}
			return b
		}
		lineA = append(lineA, parent)
		mA[parent] = struct{}{}

		child = parent
	}

	child = b
	for {
		parent, ok = t.Nodes[child]
		if !ok {
			t.lcaCache[query] = b
			return 0
		}
		if parent == child { // root
			break
		}
		if parent == a { // a is ancestor of b
			if t.cacheLCA {
				t.lcaCache[query] = a
			}
			return a
		}
		if _, ok = mA[parent]; ok {
			if t.cacheLCA {
				t.lcaCache[query] = parent
			}
			return parent
		}

		child = parent
	}
	return t.rootNode
}

func pack2uint32(a uint32, b uint32) uint64 {
	if a < b {
		return (uint64(a) << 32) | uint64(b)
	}
	return (uint64(b) << 32) | uint64(a)
}
