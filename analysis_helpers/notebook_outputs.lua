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

local function decode64(text)
  local alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/'
  local out, bits, value = {}, 0, 0
  for c in text:gmatch('[A-Za-z0-9+/]') do
    value = value * 64 + alphabet:find(c, 1, true) - 1
    bits = bits + 6
    if bits >= 8 then
      bits = bits - 8
      out[#out + 1] = string.char(math.floor(value / 2^bits))
      value = value % 2^bits
    end
  end
  return table.concat(out)
end

local function read_file(path)
  local file = io.open(path, 'r')
  if not file then return nil end
  local text = file:read('*a')
  file:close()
  return text
end

local function chunks(text)
  local all, visible = {}, {}
  for header, code in text:gmatch('```(%{r[^\n]*%})\n(.-)\n```') do
    all[#all + 1] = header .. '\n' .. code
    if not header:match('include%s*=%s*FALSE') then visible[#visible + 1] = code end
  end
  return table.concat(all, '\n'), visible
end

-- Read displayed code from the current notebook, not imported console records:
-- older RStudio versions can add fences or corrupt quotes in those records.
local function notebook_code()
  local input = PANDOC_STATE.input_files[1]
  local path = input:gsub('%.knit%.md$', '.Rmd'):gsub('%.utf8%.md$', '.Rmd')
  local text = assert(read_file(path), 'Cannot read notebook source: ' .. path)
  local signature, codes = chunks(text)
  local saved = {}
  -- A fresh RStudio session may have no cached outputs. Recover from the
  -- companion HTML only when every R chunk (including setup) is unchanged.
  local html = read_file(path:gsub('%.Rmd$', '.nb.html'))
  local encoded = html and html:match('<div id="rmd%-source%-code">(.-)</div>')
  if encoded and chunks(decode64(encoded)) == signature then
    for block in html:gmatch('<!%-%- rnb%-chunk%-begin %-%->(.-)<!%-%- rnb%-chunk%-end %-%->') do
      saved[#saved + 1] = block:match('<!%-%- rnb%-source%-end %-%->(.*)') or ''
    end
  end
  return codes, saved
end

function Pandoc(doc)
  local codes, saved = notebook_code()
  local output, stack, chunk, source = pandoc.List(), {}, 0, false
  local has_output = false
  for i, block in ipairs(doc.blocks) do
    local text = raw(block)
    local keep = true
    if text:match('^<!%-%- rnb%-chunk%-begin') then
      chunk = chunk + 1
      has_output = false
    elseif text:match('^<!%-%- rnb%-chunk%-end') then
      if not has_output and saved[chunk] and saved[chunk]:match('%S') then
        output:insert(pandoc.RawBlock('html', saved[chunk]))
      end
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
    elseif not source and not text:match('^<!%-%- rnb%-') then
      has_output = true
    end
    if keep then output:insert(block) end
  end
  doc.blocks = output
  return doc
end
