-- Author: Daniel Franks
-- Keep RStudio's saved outputs portable. Some versions nest HTML/source records
-- inside console-output records, causing duplicate results on the next import.
local function base64(text)
  local alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/'
  local out = {}
  for i = 1, #text, 3 do
    local a, b, c = text:byte(i, i + 2)
    local n = a * 65536 + (b or 0) * 256 + (c or 0)
    out[#out + 1] = alphabet:sub(math.floor(n / 262144) % 64 + 1, math.floor(n / 262144) % 64 + 1)
    out[#out + 1] = alphabet:sub(math.floor(n / 4096) % 64 + 1, math.floor(n / 4096) % 64 + 1)
    out[#out + 1] = b and alphabet:sub(math.floor(n / 64) % 64 + 1, math.floor(n / 64) % 64 + 1) or '='
    out[#out + 1] = c and alphabet:sub(n % 64 + 1, n % 64 + 1) or '='
  end
  return table.concat(out)
end

local function raw(block)
  return block and block.t == 'RawBlock' and block.format == 'html' and block.text or ''
end

-- Read displayed code from the current notebook, not imported console records:
-- older RStudio versions can add fences or corrupt quotes in those records.
local function notebook_code()
  local input = PANDOC_STATE.input_files[1]
  local path = input:gsub('%.knit%.md$', '.Rmd'):gsub('%.utf8%.md$', '.Rmd')
  local file = assert(io.open(path, 'r'), 'Cannot read notebook source: ' .. path)
  local chunks, lines, header = {}, nil, nil
  for line in file:lines() do
    if not lines and line:match('^```%{r[%s,}]') then
      header, lines = line, {}
    elseif lines and line:match('^```%s*$') then
      if not header:match('include%s*=%s*FALSE') then
        chunks[#chunks + 1] = table.concat(lines, '\n')
      end
      lines, header = nil, nil
    elseif lines then
      lines[#lines + 1] = line
    end
  end
  file:close()
  return chunks
end

function Pandoc(doc)
  local codes = notebook_code()
  local output, stack, chunk, source = pandoc.List(), {}, 0, false
  for i, block in ipairs(doc.blocks) do
    local text = raw(block)
    local keep = true
    if text:match('^<!%-%- rnb%-chunk%-begin') then
      chunk = chunk + 1
    elseif text:match('^<!%-%- rnb%-output%-begin') then
      local next_text = raw(doc.blocks[i + 1])
      local nested = next_text:match('<!%-%- rnb%-source%-begin') or
        next_text:match('<!%-%- rnb%-html%-begin') or
        next_text:match('<!%-%- rnb%-plot%-begin') or
        next_text:match('<!%-%- rnb%-frame%-begin')
      stack[#stack + 1] = not not nested
      keep = not nested
    elseif text:match('^<!%-%- rnb%-output%-end') then
      keep = not table.remove(stack)
    elseif text:match('^<!%-%- rnb%-source%-begin') then
      source = true
      local code = assert(codes[chunk], 'Notebook chunk does not match saved output')
      block = pandoc.RawBlock('html', '<!-- rnb-source-begin ' ..
        base64(pandoc.json.encode({data = code})) .. ' -->')
    elseif text:match('^<!%-%- rnb%-source%-end') then
      source = false
    elseif source and block.t == 'CodeBlock' then
      block.text = codes[chunk]
    end
    if keep then output:insert(block) end
  end
  doc.blocks = output
  return doc
end
