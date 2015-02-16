#!/usr/bin/env ruby
require 'test/unit'
require 'PDist'
require_relative 'Make_them_shorter'

class Measures < Test::Unit::TestCase
	def setup 
		@original = [1, 2, 3, 4, 5, 6]
		@perm = [5, 2, 3, 4, 1, 6]
	end
	def test_distance
		assert_equal(p, [3, 4])
	end
end 