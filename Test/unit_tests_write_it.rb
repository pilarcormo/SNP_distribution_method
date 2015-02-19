#encoding: utf-8
require_relative '../lib/write_it'
require 'test/unit'

class TestWriteIt < Test::Unit::TestCase

	def setup
		@file = "test/ratio_values.txt"
	end

	def test_file_to_array
		contents = WriteIt.file_to_array(@file)
		assert_kind_of(Array, contents)
		assert_kind_of(String, contents[0])
	end

	def test_file_to_ints_array
		contents = WriteIt.file_to_ints_array(@file)
		assert_kind_of(Array, contents)
		assert_kind_of(Integer, contents[0])
	end

	def test_file_to_floats_array
		contents = WriteIt.file_to_floats_array(@file)
		assert_kind_of(Array, contents)
		assert_kind_of(Float, contents[0])
	end

end

