#pragma once

#include <string_view>

namespace ms::file
{
	void replace_file_content(const std::string_view file_path, const std::string_view old_content, const std::string_view new_contnet);
	void remove_file(const std::string_view file_path);
} // namespace ms::file